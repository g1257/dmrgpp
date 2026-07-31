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
#include <KokkosBatched_Gemm_Decl.hpp>
#include <KokkosBatched_Gemm_TeamVector_Impl.hpp>

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
 *  Pass 2: For each ipatch:
 *     Y[ip] = BXbatch[ip] * Abatch[ip]^T
 *
 * Abatch[ip] = left-operator (xc) matrices packed column-by-column (ldA padded).
 * Bbatch[ip] = right-operator (yc) matrices packed column-by-column (ldB padded).
 * BXbatch[ip] is the intermediate work buffer (ldB padded, same col-width as Abatch).
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

    // Compact struct holding all parameters for one batched GEMM call.
    // Stored in device memory and accessed inside the kernel.
    struct GemmArgs {
        int       m, n, k;             // GEMM dimensions
        int       lda, ldb, ldc;       // leading dimensions (column-major strides)
        long long a_off, b_off, c_off; // element offsets into device flat arrays
    };

    using DevArgsView = Kokkos::View<GemmArgs*, DevMemSpace>;

    static const int ialign_ = 32;

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

        // --- Zero output and work buffers -------------------------------------------
        Kokkos::deep_copy(d_vout_, KS(0));
        Kokkos::deep_copy(d_flatBXbatch_, KS(0));

        DevExecSpace exec;

        // --- Pass 1: BXbatch[ip] = Bbatch[ip] * X[jp]  (NoT x NoT) ----------------
        {
            const DevScalView flatBbatch  = d_flatBbatch_;
            const DevScalView vin_dev     = d_vin_;
            const DevScalView flatBXbatch = d_flatBXbatch_;
            const DevArgsView args        = d_pass1_;

            using MemberType = typename Kokkos::TeamPolicy<DevExecSpace>::member_type;
            Kokkos::parallel_for(
                "BatchedGemmKokkos_Pass1",
                Kokkos::TeamPolicy<DevExecSpace>(exec, static_cast<int>(nbatch1_),
                                                 Kokkos::AUTO, Kokkos::AUTO),
                KOKKOS_LAMBDA(const MemberType& member) {
                    const int       i  = member.league_rank();
                    const GemmArgs& ag = args(i);

                    using UV = Kokkos::View<KS**, Kokkos::LayoutLeft, DevMemSpace,
                                            Kokkos::MemoryUnmanaged>;

                    // A = Bbatch[ip] block: (lda, k) padded -> subview (m, k) active rows
                    UV A_full(flatBbatch.data() + ag.a_off, ag.lda, ag.k);
                    auto A = Kokkos::subview(A_full,
                                 Kokkos::make_pair(0, ag.m), Kokkos::ALL());

                    // B = X[jp]: (k, n) exact -- no padding for X slice
                    UV B(vin_dev.data() + ag.b_off, ag.ldb, ag.n);

                    // C = BXbatch[ip] block: (ldc, n) padded -> subview (m, n) active rows
                    UV C_full(flatBXbatch.data() + ag.c_off, ag.ldc, ag.n);
                    auto C = Kokkos::subview(C_full,
                                 Kokkos::make_pair(0, ag.m), Kokkos::ALL());

                    KokkosBatched::TeamVectorGemm<
                        MemberType,
                        KokkosBatched::Trans::NoTranspose,
                        KokkosBatched::Trans::NoTranspose,
                        KokkosBatched::Algo::Gemm::Blocked>::
                        invoke(member, KS(1), A, B, KS(0), C);
                });
        }

        // --- Pass 2: Y[ip] = BXbatch[ip] * Abatch[ip]^T  (NoT x T) ---------------
        {
            const DevScalView flatBXbatch = d_flatBXbatch_;
            const DevScalView flatAbatch  = d_flatAbatch_;
            const DevScalView vout_dev    = d_vout_;
            const DevArgsView args        = d_pass2_;

            using MemberType = typename Kokkos::TeamPolicy<DevExecSpace>::member_type;
            Kokkos::parallel_for(
                "BatchedGemmKokkos_Pass2",
                Kokkos::TeamPolicy<DevExecSpace>(exec, static_cast<int>(nbatch2_),
                                                 Kokkos::AUTO, Kokkos::AUTO),
                KOKKOS_LAMBDA(const MemberType& member) {
                    const int       i  = member.league_rank();
                    const GemmArgs& ag = args(i);

                    if (ag.k == 0)
                        return; // no connections for this ipatch; Y[ip] already zero

                    using UV = Kokkos::View<KS**, Kokkos::LayoutLeft, DevMemSpace,
                                            Kokkos::MemoryUnmanaged>;

                    // A = BXbatch[ip]: (lda, k) padded -> subview (m, k)
                    UV A_full(flatBXbatch.data() + ag.a_off, ag.lda, ag.k);
                    auto A = Kokkos::subview(A_full,
                                 Kokkos::make_pair(0, ag.m), Kokkos::ALL());

                    // B = Abatch[ip]: (ldb, k) padded -> subview (n, k); will be transposed
                    UV B_full(flatAbatch.data() + ag.b_off, ag.ldb, ag.k);
                    auto B = Kokkos::subview(B_full,
                                 Kokkos::make_pair(0, ag.n), Kokkos::ALL());

                    // C = Y[ip]: (ldc, n) exact -- no padding for output
                    UV C(vout_dev.data() + ag.c_off, ag.ldc, ag.n);

                    KokkosBatched::TeamVectorGemm<
                        MemberType,
                        KokkosBatched::Trans::NoTranspose,
                        KokkosBatched::Trans::Transpose,
                        KokkosBatched::Algo::Gemm::Blocked>::
                        invoke(member, KS(1), A, B, KS(0), C);
                });
        }

        // --- D2H: copy vout back to host -------------------------------------------
        {
            using HV = Kokkos::View<KS*, Kokkos::HostSpace, Kokkos::MemoryUnmanaged>;
            HV hv(reinterpret_cast<KS*>(vout.data()), totalXY);
            Kokkos::deep_copy(exec, hv, d_vout_);
            exec.fence();
        }
    }

private:

    static int iceil(int x, int n) { return (x + n - 1) / n; }

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

        // Per-patch sizes and padded leading dimensions.
        VectorSizeType lps(npatches, 0);        // left patch size
        VectorSizeType rps(npatches, 0);        // right patch size
        VectorSizeType ldA(npatches, 0);        // padded lps -> Abatch leading dim
        VectorSizeType ldB(npatches, 0);        // padded rps -> Bbatch / BXbatch leading dim
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
            ldA[ip] = static_cast<SizeType>(ialign_ * iceil(static_cast<int>(lps[ip]), ialign_));
            ldB[ip] = static_cast<SizeType>(ialign_ * iceil(static_cast<int>(rps[ip]), ialign_));
        }

        xyStart[0] = 0;
        for (SizeType ip = 0; ip < npatches; ++ip)
            xyStart[ip + 1] = xyStart[ip] + lps[ip] * rps[ip];

        // Column widths:
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

        // Flat offsets: each patch ip occupies a contiguous slice.
        VectorSizeType AbatchOff(npatches + 1, 0);
        VectorSizeType BbatchOff(npatches + 1, 0);
        VectorSizeType BXbatchOff(npatches + 1, 0); // ldB rows x AbatchCols columns

        for (SizeType ip = 0; ip < npatches; ++ip) {
            AbatchOff[ip + 1]  = AbatchOff[ip]  + ldA[ip] * AbatchCols[ip];
            BbatchOff[ip + 1]  = BbatchOff[ip]  + ldB[ip] * BbatchCols[ip];
            BXbatchOff[ip + 1] = BXbatchOff[ip] + ldB[ip] * AbatchCols[ip];
        }

        const SizeType totalAbatch  = AbatchOff[npatches];
        const SizeType totalBbatch  = BbatchOff[npatches];
        const SizeType totalBXbatch = BXbatchOff[npatches];

        // Host packing buffers (zeroed).
        std::vector<ComplexOrRealType> h_flatAbatch(totalAbatch, ComplexOrRealType(0));
        std::vector<ComplexOrRealType> h_flatBbatch(totalBbatch, ComplexOrRealType(0));

        // Build GEMM arg lists while packing matrices.
        std::vector<GemmArgs> pass1_args;
        std::vector<GemmArgs> pass2_args;
        pass2_args.reserve(npatches);

        for (SizeType ip = 0; ip < npatches; ++ip) {
            long long colA = 0; // column cursor in Abatch[ip]  (= BXbatch[ip])
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

                    const int mA  = static_cast<int>(lps[ip]);
                    const int nAk = static_cast<int>(lps[jp]); // A cols = BX cols per op
                    const int mB  = static_cast<int>(rps[ip]);
                    const int nBk = static_cast<int>(rps[jp]); // B cols = X rows

                    // Pack A (left operator) into Abatch[ip] at column colA
                    lacpy(Aref.ptr, mA, nAk, Aref.rows,
                          h_flatAbatch.data()
                              + AbatchOff[ip] + colA * static_cast<long long>(ldA[ip]),
                          static_cast<int>(ldA[ip]));

                    // Pack B (right operator) into Bbatch[ip] at column colB
                    lacpy(Bref.ptr, mB, nBk, Bref.rows,
                          h_flatBbatch.data()
                              + BbatchOff[ip] + colB * static_cast<long long>(ldB[ip]),
                          static_cast<int>(ldB[ip]));

                    // Pass 1 GEMM: C(mB, nAk) = Bbatch_block(mB, nBk) * X_jp(nBk, nAk)
                    GemmArgs a1;
                    a1.m   = mB;
                    a1.n   = nAk;
                    a1.k   = nBk;
                    a1.lda = static_cast<int>(ldB[ip]);   // Bbatch leading dim
                    a1.ldb = nBk;                          // X[jp] no padding: ld = rps[jp]
                    a1.ldc = static_cast<int>(ldB[ip]);   // BXbatch leading dim
                    a1.a_off = static_cast<long long>(BbatchOff[ip])
                               + colB * static_cast<long long>(ldB[ip]);
                    a1.b_off = static_cast<long long>(xyStart[jp]);
                    a1.c_off = static_cast<long long>(BXbatchOff[ip])
                               + colA * static_cast<long long>(ldB[ip]);
                    pass1_args.push_back(a1);

                    colA += nAk;
                    colB += nBk;
                }
            }

            // Pass 2 GEMM: Y[ip](mB, mA) = BXbatch[ip](mB, k) * Abatch[ip](mA, k)^T
            const int totalCols = static_cast<int>(colA); // = AbatchCols[ip]
            GemmArgs  a2;
            a2.m   = static_cast<int>(rps[ip]);
            a2.n   = static_cast<int>(lps[ip]);
            a2.k   = totalCols;
            a2.lda = static_cast<int>(ldB[ip]); // BXbatch leading dim
            a2.ldb = static_cast<int>(ldA[ip]); // Abatch leading dim
            a2.ldc = static_cast<int>(rps[ip]); // Y[ip] no padding: ld = rps[ip]
            a2.a_off = static_cast<long long>(BXbatchOff[ip]);
            a2.b_off = static_cast<long long>(AbatchOff[ip]);
            a2.c_off = static_cast<long long>(xyStart[ip]);
            pass2_args.push_back(a2);
        }

        // Allocate device arrays and upload operator matrices.
        d_flatAbatch_  = DevScalView(Kokkos::view_alloc(Kokkos::WithoutInitializing, "d_flatAbatch"),  totalAbatch  ? totalAbatch  : 1);
        d_flatBbatch_  = DevScalView(Kokkos::view_alloc(Kokkos::WithoutInitializing, "d_flatBbatch"),  totalBbatch  ? totalBbatch  : 1);
        d_flatBXbatch_ = DevScalView(Kokkos::view_alloc(Kokkos::WithoutInitializing, "d_flatBXbatch"), totalBXbatch ? totalBXbatch : 1);
        d_vin_         = DevScalView(Kokkos::view_alloc(Kokkos::WithoutInitializing, "d_vin"),  xyStart[npatches] ? xyStart[npatches] : 1);
        d_vout_        = DevScalView(Kokkos::view_alloc(Kokkos::WithoutInitializing, "d_vout"), xyStart[npatches] ? xyStart[npatches] : 1);

        nbatch1_ = pass1_args.size();
        nbatch2_ = pass2_args.size(); // == npatches

        {
            auto hA = Kokkos::create_mirror_view(Kokkos::view_alloc(Kokkos::WithoutInitializing), d_flatAbatch_);
            auto hB = Kokkos::create_mirror_view(Kokkos::view_alloc(Kokkos::WithoutInitializing), d_flatBbatch_);
            for (SizeType i = 0; i < totalAbatch; ++i)
                hA(i) = *reinterpret_cast<const KS*>(&h_flatAbatch[i]);
            for (SizeType i = 0; i < totalBbatch; ++i)
                hB(i) = *reinterpret_cast<const KS*>(&h_flatBbatch[i]);
            Kokkos::deep_copy(d_flatAbatch_, hA);
            Kokkos::deep_copy(d_flatBbatch_, hB);
        }

        d_pass1_ = DevArgsView("d_pass1", nbatch1_ ? nbatch1_ : 1);
        d_pass2_ = DevArgsView("d_pass2", nbatch2_ ? nbatch2_ : 1);
        {
            auto h1 = Kokkos::create_mirror_view(d_pass1_);
            auto h2 = Kokkos::create_mirror_view(d_pass2_);
            for (SizeType i = 0; i < nbatch1_; ++i) h1(i) = pass1_args[i];
            for (SizeType i = 0; i < nbatch2_; ++i) h2(i) = pass2_args[i];
            Kokkos::deep_copy(d_pass1_, h1);
            Kokkos::deep_copy(d_pass2_, h2);
        }

        {
            PsimagLite::OstringStream msg(std::cout.precision());
            msg() << "setup done: npatches=" << npatches
                  << " noperator=" << noperator
                  << " pass1_batches=" << nbatch1_
                  << " pass2_batches=" << nbatch2_
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

    SizeType nbatch1_ = 0; // number of pass-1 GEMMs
    SizeType nbatch2_ = 0; // number of pass-2 GEMMs (== npatches)

    DevScalView          d_flatAbatch_;   // left-operator matrices, device, persistent
    DevScalView          d_flatBbatch_;   // right-operator matrices, device, persistent
    mutable DevScalView  d_flatBXbatch_;  // intermediate BX work buffer, per-call
    mutable DevScalView  d_vin_;          // input vector on device, per-call
    mutable DevScalView  d_vout_;         // output vector on device, per-call
    DevArgsView          d_pass1_;        // pass-1 GEMM parameters, persistent
    DevArgsView          d_pass2_;        // pass-2 GEMM parameters, persistent
};

} // namespace Dmrg
#endif // BATCHED_GEMM_KOKKOS_H
