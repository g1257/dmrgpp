#ifndef BATCHED_GEMM_KOKKOS_H
#define BATCHED_GEMM_KOKKOS_H
// Don't include this file directly; use BatchedGemmInclude.hh

#include <PsimagLite/CrsMatrix.h>
#include <PsimagLite/KokkosType.h>
#include <PsimagLite/Matrix.h>
#include <PsimagLite/ProgressIndicator.h>
#include <PsimagLite/PsimagLite.h>
#include <PsimagLite/Vector.h>

#include <Kokkos_Core.hpp>
#include <Kokkos_Profiling_ScopedRegion.hpp>
#include <KokkosBlas3_gemm.hpp>

namespace Dmrg {

/*
 * Batched GEMM using KokkosKernels with a sparse-batch layout.
 *
 * Algorithm (mirrors apply_Htarget_sparse_not_ready.cpp):
 *
 *  Pass 1: For each non-zero triple (ipatch, jpatch, ioperator):
 *     BXbatch[ip][:, colBX:colBX+lps[jp]] =
 *         Bbatch[ip][:, colB:colB+rps[jp]] * X[jp]
 *
 *  Pass 2: For each ipatch (with connections):
 *     Y[ip] = BXbatch[ip] * Abatch[ip]^T
 *
 * Abatch[ip] = left-operator (xc) matrices packed column-by-column (no padding).
 * Bbatch[ip] = right-operator (yc) matrices packed column-by-column (no padding).
 * BXbatch[ip] is the intermediate work buffer (no padding, same col-width as Abatch).
 *
 * No padding is used so that unmanaged Kokkos::LayoutLeft views can be created
 * directly from device pointers, enabling dispatch to rocBLAS/cuBLAS via
 * KokkosBlas::gemm.
 *
 * All Abatch and Bbatch data live on the device for the lifetime of the object.
 * BXbatch, vin, vout are (re-)written on each matrixVector call.
 */

template <typename InitKronType>
class BatchedGemmKokkos {

    using ArrayOfMatStructType    = typename InitKronType::ArrayOfMatStructType;
    using GenIjPatchType          = typename InitKronType::GenIjPatchType;
    using MatrixDenseOrSparseType = typename ArrayOfMatStructType::MatrixDenseOrSparseType;
    using VectorType              = typename MatrixDenseOrSparseType::VectorType;
    using ComplexOrRealType       = typename VectorType::value_type;
    using MatrixType              = PsimagLite::Matrix<ComplexOrRealType>;
    using VectorMatrixType        = typename PsimagLite::Vector<MatrixType*>::Type;
    using VectorSizeType          = PsimagLite::Vector<SizeType>::Type;

    // Kokkos types
    using KS           = typename PsimagLite::KokkosType<ComplexOrRealType>::type;
    using DevExecSpace = Kokkos::DefaultExecutionSpace;
    using DevMemSpace  = typename DevExecSpace::memory_space;
    using DevScalView  = Kokkos::View<KS*, DevMemSpace>;

    // Parameters for one KokkosBlas::gemm call.
    // With no padding: leading dims = m (for A,C) or k (for B in pass1) or n (for B in pass2).
    struct GemmArgs {
        int       m, n, k;             // GEMM dimensions
        long long a_off, b_off, c_off; // element offsets into device flat arrays
    };

public:

    BatchedGemmKokkos(const InitKronType& initKron)
        : initKron_(initKron)
        , progress_("BatchedGemmKokkos")
    {
        if (!enabled())
            return;
        setup_();
    }

    ~BatchedGemmKokkos()
    {
        for (SizeType i = 0; i < garbage_.size(); ++i) {
            delete garbage_[i];
            garbage_[i] = nullptr;
        }
    }

    bool enabled() const { return initKron_.params().options.isSet("BatchedGemm"); }

    void matrixVector(VectorType& vout, const VectorType& vin) const
    {
        assert(enabled());

        Kokkos::Profiling::ScopedRegion region("BatchedGemmKokkos::matrixVector");

        const SizeType totalXY = vin.size();
        assert(vout.size() == totalXY);
        assert(d_vin_.extent(0) == totalXY);

        // --- H2D: copy vin to device -----------------------------------------------
        {
            using HV = Kokkos::View<const KS*, Kokkos::HostSpace, Kokkos::MemoryUnmanaged>;
            HV hv(reinterpret_cast<const KS*>(vin.data()), totalXY);
            Kokkos::deep_copy(d_vin_, hv);
        }

        // Zero d_vout_ so patches with no connections produce zero output.
        // d_flatBXbatch_ is NOT zeroed: pass1 uses beta=0 (overwrites all used blocks).
        Kokkos::deep_copy(d_vout_, KS(0));

        DevExecSpace exec;

        using UV = Kokkos::View<KS**, Kokkos::LayoutLeft, DevMemSpace, Kokkos::MemoryUnmanaged>;

        // --- Pass 1: BXbatch[ip] block = Bbatch[ip] block * X[jp]  (N x N) ---------
        // Each call: C(m, n) = A(m, k) * B(k, n), beta=0 overwrites C.
        //   A = Bbatch block  (m=rps[ip], k=rps[jp]),  ld = rps[ip]
        //   B = X[jp]         (k=rps[jp], n=lps[jp]),  ld = rps[jp]
        //   C = BXbatch block (m=rps[ip], n=lps[jp]),  ld = rps[ip]
        for (const GemmArgs& ag : h_pass1_) {
            UV A(d_flatBbatch_.data()  + ag.a_off, ag.m, ag.k);
            UV B(d_vin_.data()         + ag.b_off, ag.k, ag.n);
            UV C(d_flatBXbatch_.data() + ag.c_off, ag.m, ag.n);
            KokkosBlas::gemm(exec, "N", "N", KS(1), A, B, KS(0), C);
        }

        exec.fence();

        // --- Pass 2: Y[ip] = BXbatch[ip] * Abatch[ip]^T  (N x T) ------------------
        // Each call: C(m, n) = A(m, k) * B(n, k)^T, beta=0 overwrites C.
        //   A = BXbatch[ip]  (m=rps[ip], k=AbatchCols[ip]), ld = rps[ip]
        //   B = Abatch[ip]   (n=lps[ip], k=AbatchCols[ip]), ld = lps[ip]  (transposed)
        //   C = Y[ip]        (m=rps[ip], n=lps[ip]),         ld = rps[ip]
        for (const GemmArgs& ag : h_pass2_) {
            if (ag.k == 0)
                continue; // no connections for this ipatch; Y[ip] stays zero
            UV A(d_flatBXbatch_.data() + ag.a_off, ag.m, ag.k);
            UV B(d_flatAbatch_.data()  + ag.b_off, ag.n, ag.k);
            UV C(d_vout_.data()        + ag.c_off, ag.m, ag.n);
            KokkosBlas::gemm(exec, "N", "T", KS(1), A, B, KS(0), C);
        }

        exec.fence();

        // --- D2H: copy vout back to host -------------------------------------------
        {
            using HV = Kokkos::View<KS*, Kokkos::HostSpace, Kokkos::MemoryUnmanaged>;
            HV hv(reinterpret_cast<KS*>(vout.data()), totalXY);
            Kokkos::deep_copy(hv, d_vout_);
        }
    }

private:

    // Dense matrix reference: pointer + source dimensions (no padding).
    // For sparse matrices, the expansion is heap-allocated and owned by garbage_.
    struct MatRef {
        ComplexOrRealType* ptr;
        int                rows; // = source leading dimension
        int                cols;
    };

    MatRef getMatRef(const MatrixDenseOrSparseType& mat)
    {
        if (mat.isZero())
            return { nullptr, 0, 0 };
        if (!mat.isDense()) {
            MatrixType* dense = new MatrixType();
            crsMatrixToFullMatrix(*dense, mat.sparse());
            garbage_.push_back(dense);
            return { &(*dense)(0, 0),
                     static_cast<int>(dense->rows()),
                     static_cast<int>(dense->cols()) };
        }
        const MatrixType& d = mat.dense();
        return { const_cast<ComplexOrRealType*>(&d(0, 0)),
                 static_cast<int>(d.rows()),
                 static_cast<int>(d.cols()) };
    }

    // Column-major copy: src(m, n) with lds -> dst(m, n) with ldd.
    static void lacpy(const ComplexOrRealType* src, int m, int n, int lds,
                      ComplexOrRealType* dst, int ldd)
    {
        for (int j = 0; j < n; ++j)
            for (int i = 0; i < m; ++i)
                dst[i + static_cast<long long>(j) * ldd]
                    = src[i + static_cast<long long>(j) * lds];
    }

    void setup_()
    {
        Kokkos::Profiling::ScopedRegion region("BatchedGemmKokkos::setup");

        const SizeType npatches  = initKron_.numberOfPatches(InitKronType::OLD);
        const SizeType noperator = initKron_.connections();

        // Per-patch sizes — NO padding: leading dim == patch size.
        VectorSizeType lps(npatches, 0);         // left patch size
        VectorSizeType rps(npatches, 0);         // right patch size
        VectorSizeType xyStart(npatches + 1, 0); // patch start offset in vin/vout

        for (SizeType ip = 0; ip < npatches; ++ip) {
            const SizeType lg = initKron_.patch(InitKronType::NEW, GenIjPatchType::LEFT)[ip];
            const SizeType rg = initKron_.patch(InitKronType::NEW, GenIjPatchType::RIGHT)[ip];
            const int L1 = initKron_.lrs(InitKronType::NEW).left().partition(lg);
            const int L2 = initKron_.lrs(InitKronType::NEW).left().partition(lg + 1);
            const int R1 = initKron_.lrs(InitKronType::NEW).right().partition(rg);
            const int R2 = initKron_.lrs(InitKronType::NEW).right().partition(rg + 1);
            lps[ip] = static_cast<SizeType>(L2 - L1);
            rps[ip] = static_cast<SizeType>(R2 - R1);
        }

        xyStart[0] = 0;
        for (SizeType ip = 0; ip < npatches; ++ip)
            xyStart[ip + 1] = xyStart[ip] + lps[ip] * rps[ip];

        // Column widths for each patch's batch block.
        //   AbatchCols[ip] = sum of lps[jp] over all non-zero (jp, k) connections for ip
        //   BbatchCols[ip] = sum of rps[jp] over all non-zero (jp, k) connections for ip
        VectorSizeType AbatchCols(npatches, 0);
        VectorSizeType BbatchCols(npatches, 0);

        for (SizeType ip = 0; ip < npatches; ++ip) {
            for (SizeType jp = 0; jp < npatches; ++jp) {
                for (SizeType k = 0; k < noperator; ++k) {
                    const MatrixDenseOrSparseType* Ap = initKron_.xc(k)(ip, jp);
                    const MatrixDenseOrSparseType* Bp = initKron_.yc(k)(ip, jp);
                    if (!Ap || !Bp) continue;
                    if (Ap->isZero() || Bp->isZero()) continue;
                    AbatchCols[ip] += lps[jp];
                    BbatchCols[ip] += rps[jp];
                }
            }
        }

        // Flat offsets: no padding, ld[ip] = rps[ip] or lps[ip].
        //   Abatch[ip]:  lps[ip] rows x AbatchCols[ip] cols
        //   Bbatch[ip]:  rps[ip] rows x BbatchCols[ip] cols
        //   BXbatch[ip]: rps[ip] rows x AbatchCols[ip] cols
        VectorSizeType AbatchOff(npatches + 1, 0);
        VectorSizeType BbatchOff(npatches + 1, 0);
        VectorSizeType BXbatchOff(npatches + 1, 0);

        for (SizeType ip = 0; ip < npatches; ++ip) {
            AbatchOff[ip + 1]  = AbatchOff[ip]  + lps[ip] * AbatchCols[ip];
            BbatchOff[ip + 1]  = BbatchOff[ip]  + rps[ip] * BbatchCols[ip];
            BXbatchOff[ip + 1] = BXbatchOff[ip] + rps[ip] * AbatchCols[ip];
        }

        const SizeType totalAbatch  = AbatchOff[npatches];
        const SizeType totalBbatch  = BbatchOff[npatches];
        const SizeType totalBXbatch = BXbatchOff[npatches];

        // Host packing buffers (zeroed).
        std::vector<ComplexOrRealType> h_flatAbatch(totalAbatch, ComplexOrRealType(0));
        std::vector<ComplexOrRealType> h_flatBbatch(totalBbatch, ComplexOrRealType(0));

        // Build host GEMM arg lists while packing operator matrices.
        h_pass1_.clear();
        h_pass2_.clear();
        h_pass2_.reserve(npatches);

        for (SizeType ip = 0; ip < npatches; ++ip) {
            long long colA = 0; // column cursor in Abatch[ip] and BXbatch[ip]
            long long colB = 0; // column cursor in Bbatch[ip]

            for (SizeType jp = 0; jp < npatches; ++jp) {
                for (SizeType k = 0; k < noperator; ++k) {
                    const MatrixDenseOrSparseType* Amat = initKron_.xc(k)(ip, jp);
                    const MatrixDenseOrSparseType* Bmat = initKron_.yc(k)(ip, jp);
                    if (!Amat || !Bmat) continue;
                    if (Amat->isZero() || Bmat->isZero()) continue;

                    MatRef Aref = getMatRef(*Amat); // left operator (xc)
                    MatRef Bref = getMatRef(*Bmat); // right operator (yc)
                    assert(Aref.ptr && Bref.ptr);

                    const int mA  = static_cast<int>(lps[ip]); // Abatch rows
                    const int nAk = static_cast<int>(lps[jp]); // Abatch cols per op
                    const int mB  = static_cast<int>(rps[ip]); // Bbatch / BXbatch rows
                    const int nBk = static_cast<int>(rps[jp]); // Bbatch cols per op

                    // Pack A (left operator) into Abatch[ip] at column colA.
                    // Abatch[ip]: lps[ip] rows, ld = lps[ip] (no padding).
                    lacpy(Aref.ptr, mA, nAk, Aref.rows,
                          h_flatAbatch.data()
                              + AbatchOff[ip] + colA * static_cast<long long>(lps[ip]),
                          mA);

                    // Pack B (right operator) into Bbatch[ip] at column colB.
                    // Bbatch[ip]: rps[ip] rows, ld = rps[ip] (no padding).
                    lacpy(Bref.ptr, mB, nBk, Bref.rows,
                          h_flatBbatch.data()
                              + BbatchOff[ip] + colB * static_cast<long long>(rps[ip]),
                          mB);

                    // Pass 1 GEMM: C(mB, nAk) = A(mB, nBk) * B(nBk, nAk)
                    //   A = Bbatch block, ld = rps[ip] = mB
                    //   B = X[jp],        ld = rps[jp] = nBk
                    //   C = BXbatch block, ld = rps[ip] = mB
                    GemmArgs a1;
                    a1.m     = mB;
                    a1.n     = nAk;
                    a1.k     = nBk;
                    a1.a_off = static_cast<long long>(BbatchOff[ip])
                               + colB * static_cast<long long>(rps[ip]);
                    a1.b_off = static_cast<long long>(xyStart[jp]);
                    a1.c_off = static_cast<long long>(BXbatchOff[ip])
                               + colA * static_cast<long long>(rps[ip]);
                    h_pass1_.push_back(a1);

                    colA += nAk;
                    colB += nBk;
                }
            }

            // Pass 2 GEMM: Y[ip](mB, mA) = BXbatch[ip](mB, k) * Abatch[ip](mA, k)^T
            //   A = BXbatch[ip], ld = rps[ip] = mB
            //   B = Abatch[ip],  ld = lps[ip] = mA  (transposed)
            //   C = Y[ip],       ld = rps[ip] = mB
            const int totalCols = static_cast<int>(colA); // = AbatchCols[ip]
            GemmArgs  a2;
            a2.m     = static_cast<int>(rps[ip]);
            a2.n     = static_cast<int>(lps[ip]);
            a2.k     = totalCols;
            a2.a_off = static_cast<long long>(BXbatchOff[ip]);
            a2.b_off = static_cast<long long>(AbatchOff[ip]);
            a2.c_off = static_cast<long long>(xyStart[ip]);
            h_pass2_.push_back(a2);
        }

        // Allocate device arrays and upload operator matrices.
        d_flatAbatch_  = DevScalView("d_flatAbatch",  totalAbatch  ? totalAbatch  : 1);
        d_flatBbatch_  = DevScalView("d_flatBbatch",  totalBbatch  ? totalBbatch  : 1);
        d_flatBXbatch_ = DevScalView("d_flatBXbatch", totalBXbatch ? totalBXbatch : 1);
        d_vin_         = DevScalView("d_vin",  xyStart[npatches] ? xyStart[npatches] : 1);
        d_vout_        = DevScalView("d_vout", xyStart[npatches] ? xyStart[npatches] : 1);

        {
            auto hA = Kokkos::create_mirror_view(d_flatAbatch_);
            auto hB = Kokkos::create_mirror_view(d_flatBbatch_);
            for (SizeType i = 0; i < totalAbatch; ++i)
                hA(i) = *reinterpret_cast<const KS*>(&h_flatAbatch[i]);
            for (SizeType i = 0; i < totalBbatch; ++i)
                hB(i) = *reinterpret_cast<const KS*>(&h_flatBbatch[i]);
            Kokkos::deep_copy(d_flatAbatch_, hA);
            Kokkos::deep_copy(d_flatBbatch_, hB);
        }

        {
            PsimagLite::OstringStream msg(std::cout.precision());
            msg() << "setup done: npatches=" << npatches
                  << " noperator=" << noperator
                  << " pass1_batches=" << h_pass1_.size()
                  << " pass2_batches=" << h_pass2_.size()
                  << " Abatch=" << totalAbatch << "elems"
                  << " Bbatch=" << totalBbatch << "elems"
                  << " BXbatch=" << totalBXbatch << "elems";
            progress_.printline(msg, std::cout);
        }
    }

    // -----------------------------------------------------------------------

    const InitKronType&           initKron_;
    PsimagLite::ProgressIndicator progress_;
    mutable VectorMatrixType      garbage_; // owns sparse->dense expansions (freed in dtor)

    std::vector<GemmArgs>        h_pass1_; // pass-1 GEMM parameters (host, persistent)
    std::vector<GemmArgs>        h_pass2_; // pass-2 GEMM parameters (host, persistent)

    DevScalView          d_flatAbatch_;   // left-operator matrices, device, persistent
    DevScalView          d_flatBbatch_;   // right-operator matrices, device, persistent
    mutable DevScalView  d_flatBXbatch_;  // intermediate BX work buffer, per-call
    mutable DevScalView  d_vin_;          // input vector on device, per-call
    mutable DevScalView  d_vout_;         // output vector on device, per-call
};

} // namespace Dmrg
#endif // BATCHED_GEMM_KOKKOS_H
