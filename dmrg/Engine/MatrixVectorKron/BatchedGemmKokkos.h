#ifndef BATCHED_GEMM_KOKKOS_H
#define BATCHED_GEMM_KOKKOS_H
// Don't include this file directly; use BatchedGemmInclude.hh

#include <PsimagLite/CrsMatrix.h>
#include <PsimagLite/KokkosType.h>
#include <PsimagLite/Matrix.h>
#include <PsimagLite/ProgressIndicator.h>
#include <PsimagLite/PsimagLite.h>
#include <PsimagLite/Vector.h>

#include <KokkosBatched_Gemm_Decl.hpp>
#include <KokkosBatched_Gemm_TeamVector_Impl.hpp>
#include <Kokkos_Core.hpp>
#include <Kokkos_Profiling_ScopedRegion.hpp>

#include <algorithm>
#include <memory>

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

template <typename InitKronType> class BatchedGemmKokkos {

	using ArrayOfMatStructType    = typename InitKronType::ArrayOfMatStructType;
	using GenIjPatchType          = typename InitKronType::GenIjPatchType;
	using MatrixDenseOrSparseType = typename ArrayOfMatStructType::MatrixDenseOrSparseType;
	using VectorType              = typename MatrixDenseOrSparseType::VectorType;
	using ComplexOrRealType       = typename VectorType::value_type;
	using MatrixType              = PsimagLite::Matrix<ComplexOrRealType>;
	using VectorMatrixType = typename PsimagLite::Vector<std::unique_ptr<MatrixType>>::Type;
	using VectorSizeType   = PsimagLite::Vector<SizeType>::Type;

	// Kokkos types
	using KokkosScalar   = typename PsimagLite::KokkosType<ComplexOrRealType>::type;
	using ExecutionSpace = Kokkos::DefaultExecutionSpace;
	using MemorySpace    = typename ExecutionSpace::memory_space;
	using ScalarView     = Kokkos::View<KokkosScalar*, MemorySpace>;

	// Needs to be public so we can use lambdas with Cuda.
public:

	// Compact struct holding all parameters for one batched GEMM call.
	// Stored in device memory and accessed inside the kernel.
	struct GemmArgs {
		int       m, n, k; // GEMM dimensions
		int       lda, ldb, ldc; // leading dimensions (column-major strides)
		long long a_off, b_off, c_off; // element offsets into device flat arrays
	};

private:

	using DevArgsView = Kokkos::View<GemmArgs*, MemorySpace>;

	static const int ialign_ = 32;

	// Pass 2 GEMMs have very few batches (one per patch, == npatches, which
	// can be as low as single digits) but each has a very long contraction
	// (k) dimension, so a single GEMM per patch leaves the GPU with far too
	// few concurrent thread-blocks to reach good occupancy. We split the k
	// dimension of each patch's GEMM into independent slices and accumulate
	// their partial results directly into d_vout_ via atomic adds (see
	// matrixVector()). This multiplies the number of pass-2 batches (and
	// thus concurrent thread-blocks) without changing the per-GEMM output
	// shape (so the favorable Blocked-algo memory access pattern is
	// preserved).
	//
	// Splitting into a *fixed count* of chunks per patch (tried first)
	// under-performs badly when patches have very different k (contraction)
	// sizes: every patch still gets the same number of chunks, so patches
	// with large k get large (slow) chunks while patches with small k get
	// tiny (fast) chunks. Since all chunks across all patches are launched
	// in the same kernel, the wall-clock time is dominated by the few
	// largest chunks while most thread-blocks finish almost instantly and
	// leave the GPU idle (observed via ncu: waves/SM as low as 0.12-0.30,
	// achieved occupancy 5-8%, despite ~56% theoretical occupancy).
	//
	// Instead we split by a *fixed column width* (kPass2ColChunk_): each
	// patch gets ceil(k_ip / kPass2ColChunk_) chunks, so patches with larger
	// k automatically get proportionally more (similarly-sized) chunks.
	// This balances per-team work across the whole grid regardless of how
	// unevenly k is distributed across patches.
	//
	// Reduction strategy: an earlier version of this split-k scheme wrote
	// each chunk's partial (m x n) result into its own private slice of a
	// large global scratch buffer (sized chunks_per_patch * m * n, summed
	// over all patches), then ran a separate reduction kernel to sum the
	// slices into d_vout_. That scratch buffer's size grows with the
	// *product* of (number of chunks) and (patch output size), both of
	// which grow with system size -- for large problems (many kept states,
	// many sites) this blew up to 100+ GiB and caused out-of-memory
	// failures. Instead, each chunk-team now computes its partial result
	// into a small *per-team* scratch buffer backed by Kokkos' global-
	// memory scratch pool (level 1), whose total size is bounded by the
	// number of *concurrently resident* teams (a small, GPU-occupancy-
	// bounded constant) times the per-team scratch size -- NOT by the
	// total number of chunks across the whole batch. The team then
	// atomically adds its scratch buffer into d_vout_ (zeroed before
	// Pass 2). This eliminates the chunk-count-dependent memory blowup.
	//
	// However, a single patch's own output (m x n = rps[ip] * lps[ip]) can
	// itself be very large for big systems (observed: a single patch's
	// scratch requirement reached tens of GiB when multiplied across the
	// GPU's many concurrently-resident teams). To keep the per-team
	// scratch size bounded by a small constant regardless of problem size,
	// we additionally tile each patch's (m x n) output into at most
	// kPass2TileDim_ x kPass2TileDim_ blocks (in addition to the k-chunking
	// above); output tiles/rows are independent (no reduction needed
	// across tiles, only across k-chunks of the same tile), so this is a
	// straightforward extra split, analogous to Pass 1's row-splitting.
	static const int kPass2ColChunk_ = 16;
	static const int kPass2TileDim_  = 128;

	// Pass 1 batches one GEMM per non-zero (ip, jp, k) connection triple, but
	// the number of such triples (tens to low thousands, see setup_) is
	// often far smaller than what's needed to fill a modern GPU with
	// concurrent thread-blocks, even though many of those GEMMs have a
	// sizeable number of output rows (m = rps[ip], which can range into the
	// hundreds). Since output rows are fully independent (no reduction
	// needed, unlike the k-split used for Pass 2), we split each triple's m
	// dimension into row-blocks of at most kPass1RowChunk_ rows and emit one
	// GEMM per row-block, multiplying the number of pass-1 teams for
	// large-m triples while leaving small-m triples (m <= kPass1RowChunk_)
	// unsplit.
	static const int kPass1RowChunk_ = 16;

public:

	BatchedGemmKokkos(const InitKronType& initKron)
	    : initKron_(initKron)
	    , progress_("BatchedGemmKokkos")
	{
		if (!enabled())
			return;
		setup_();
	}

	bool enabled() const { return initKron_.params().options.isSet("BatchedGemm"); }

	void matrixVector(VectorType& vout, const VectorType& vin) const
	{
		assert(enabled());

		Kokkos::Profiling::ScopedRegion region("BatchedGemmKokkos::matrixVector");

		const SizeType totalXY = vin.size();
		assert(vout.size() == totalXY);
		assert(d_vin_.extent(0) == totalXY);

		ExecutionSpace exec;

		// --- H2D: copy vin to device -----------------------------------------------
		{
			using HV = Kokkos::
			    View<const KokkosScalar*, Kokkos::HostSpace, Kokkos::MemoryUnmanaged>;
			HV hv(reinterpret_cast<const KokkosScalar*>(vin.data()), totalXY);
			Kokkos::deep_copy(exec, d_vin_, hv);
		}

		// Note: d_flatBXbatch_ is NOT explicitly zeroed here. Every element
		// is written with beta=0 by TeamVectorGemmInternal inside Pass 1
		// (which zero-fills C whenever beta==0, even for k==0 slices), so an
		// upfront deep_copy zeroing would be redundant extra global-memory
		// traffic. d_vout_, however, IS explicitly zeroed just before Pass 2
		// below, since Pass 2's split-k chunks now accumulate into it via
		// atomic adds (see kPass2ColChunk_ comment above) rather than each
		// writing a disjoint region.

		// --- Pass 1: BXbatch[ip] = Bbatch[ip] * X[jp]  (NoT x NoT) ----------------
		{
			const ScalarView  flatBbatch  = d_flatBbatch_;
			const ScalarView  vin_dev     = d_vin_;
			const ScalarView  flatBXbatch = d_flatBXbatch_;
			const DevArgsView args        = d_pass1_;

			using MemberType = typename Kokkos::TeamPolicy<ExecutionSpace>::member_type;
			Kokkos::parallel_for(
			    "BatchedGemmKokkos_Pass1",
			    Kokkos::TeamPolicy<ExecutionSpace>(
			        exec, static_cast<int>(nbatch1_), Kokkos::AUTO, 8),
			    KOKKOS_LAMBDA(const MemberType& member) {
				    const int       i  = member.league_rank();
				    const GemmArgs& ag = args(i);

				    using UV = Kokkos::View<KokkosScalar**,
				                            Kokkos::LayoutLeft,
				                            MemorySpace,
				                            Kokkos::MemoryUnmanaged>;

				    // A = Bbatch[ip] block: (lda, k) padded -> subview (m, k)
				    // active rows
				    UV   A_full(flatBbatch.data() + ag.a_off, ag.lda, ag.k);
				    auto A = Kokkos::subview(
				        A_full, Kokkos::make_pair(0, ag.m), Kokkos::ALL());

				    // B = X[jp]: (k, n) exact -- no padding for X slice
				    UV B(vin_dev.data() + ag.b_off, ag.ldb, ag.n);

				    // C = BXbatch[ip] block: (ldc, n) padded -> subview (m, n)
				    // active rows
				    UV   C_full(flatBXbatch.data() + ag.c_off, ag.ldc, ag.n);
				    auto C = Kokkos::subview(
				        C_full, Kokkos::make_pair(0, ag.m), Kokkos::ALL());

				    KokkosBatched::TeamVectorGemm<
				        MemberType,
				        KokkosBatched::Trans::NoTranspose,
				        KokkosBatched::Trans::NoTranspose,
				        KokkosBatched::Algo::Gemm::Blocked>::invoke(member,
				                                                    KokkosScalar(1),
				                                                    A,
				                                                    B,
				                                                    KokkosScalar(0),
				                                                    C);
			    });
		}

		// --- Pass 2: Y[ip] += BXbatch[ip]_chunk * Abatch[ip]_chunk^T --------------
		// (split-k: each patch's contraction dim is divided into
		// independent column-width-bounded slices for better GPU occupancy;
		// see kPass2ColChunk_ comment above. Each chunk-team computes its
		// partial (m x n) result into a small per-team scratch buffer -- NOT
		// a global buffer sized by total chunk count -- then atomically
		// adds it into d_vout_, so d_vout_ must be zeroed first.)
		{
			Kokkos::deep_copy(exec, d_vout_, KokkosScalar(0));

			const ScalarView  flatBXbatch = d_flatBXbatch_;
			const ScalarView  flatAbatch  = d_flatAbatch_;
			const ScalarView  vout_dev    = d_vout_;
			const DevArgsView args        = d_pass2_;
			// Fixed, compile-time-bounded per-team scratch size: every
			// Pass 2 GEMM's output tile is at most kPass2TileDim_ x
			// kPass2TileDim_ (see setup_()), so this does NOT grow with
			// problem size, unlike an earlier version that sized the
			// scratch by the (unbounded) largest single patch's full
			// output.
			constexpr size_t scratchBytesPerTeam
			    = static_cast<size_t>(kPass2TileDim_) * kPass2TileDim_ * sizeof(KokkosScalar);

			using MemberType = typename Kokkos::TeamPolicy<ExecutionSpace>::member_type;
			using ScratchView = Kokkos::View<KokkosScalar**,
			                                 Kokkos::LayoutLeft,
			                                 typename ExecutionSpace::scratch_memory_space,
			                                 Kokkos::MemoryUnmanaged>;
			Kokkos::parallel_for(
			    "BatchedGemmKokkos_Pass2",
			    Kokkos::TeamPolicy<ExecutionSpace>(
			        exec, static_cast<int>(nbatch2_), Kokkos::AUTO, 8)
			        .set_scratch_size(1, Kokkos::PerTeam(static_cast<int>(scratchBytesPerTeam))),
			    KOKKOS_LAMBDA(const MemberType& member) {
				    const int       i  = member.league_rank();
				    const GemmArgs& ag = args(i);

				    // Note: we deliberately do NOT special-case ag.k == 0 here.
				    // TeamVectorGemmInternal zero-fills C whenever beta == 0,
				    // even if k == 0, so every element of this chunk's output
				    // scratch is always written by the invoke() call below.

				    using UV = Kokkos::View<KokkosScalar**,
				                            Kokkos::LayoutLeft,
				                            MemorySpace,
				                            Kokkos::MemoryUnmanaged>;

				    // A = BXbatch[ip] chunk: (lda, k) padded -> subview (m, k)
				    UV   A_full(flatBXbatch.data() + ag.a_off, ag.lda, ag.k);
				    auto A = Kokkos::subview(
				        A_full, Kokkos::make_pair(0, ag.m), Kokkos::ALL());

				    // B = Abatch[ip] chunk: (ldb, k) padded -> subview (n, k); will
				    // be transposed
				    UV   B_full(flatAbatch.data() + ag.b_off, ag.ldb, ag.k);
				    auto B = Kokkos::subview(
				        B_full, Kokkos::make_pair(0, ag.n), Kokkos::ALL());

				    // C = this chunk's partial Y[ip]: small per-team scratch
				    // buffer (m x n, contiguous), NOT a slice of a global
				    // buffer sized by total chunk count -- see comment above.
				    ScratchView C(member.team_scratch(1), ag.m, ag.n);

				    // Pass 2 GEMMs have a very "thin" output (m, n are the
				    // per-patch sizes) with a very long contraction
				    // dimension k. The Blocked algorithm reuses loaded A/B
				    // tiles across the register-blocked inner loop, which is
				    // far more memory-bandwidth efficient for this thin/long
				    // GEMM shape than assigning one thread per output
				    // element (tried Unblocked: ~3x slower here because
				    // every thread re-reads all of A/B from global memory).
				    KokkosBatched::TeamVectorGemm<
				        MemberType,
				        KokkosBatched::Trans::NoTranspose,
				        KokkosBatched::Trans::Transpose,
				        KokkosBatched::Algo::Gemm::Blocked>::invoke(member,
				                                                    KokkosScalar(1),
				                                                    A,
				                                                    B,
				                                                    KokkosScalar(0),
				                                                    C);

				    // Accumulate this chunk's partial result into the
				    // shared global output via atomic adds (multiple chunks
				    // of the same patch race on the same d_vout_ region).
				    // Note: ag.ldc is the FULL patch's row count (rps[ip]),
				    // used as the addressing stride into d_vout_ -- NOT the
				    // same as ag.m, which is only this tile's row count
				    // (see kPass2TileDim_ comment above).
				    member.team_barrier();
				    Kokkos::parallel_for(
				        Kokkos::TeamThreadRange(member, ag.n), [&](const int j) {
					        Kokkos::parallel_for(
					            Kokkos::ThreadVectorRange(member, ag.m),
					            [&](const int r) {
						            Kokkos::atomic_add(
						                &vout_dev(ag.c_off + r
						                    + static_cast<long long>(j) * ag.ldc),
						                C(r, j));
					            });
				        });
			    });
		}

		// --- D2H: copy vout back to host -------------------------------------------
		{
			using HV = Kokkos::
			    View<KokkosScalar*, Kokkos::HostSpace, Kokkos::MemoryUnmanaged>;
			HV hv(reinterpret_cast<KokkosScalar*>(vout.data()), totalXY);
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
			garbage_.emplace_back(dense);
			return { &(*dense)(0, 0),
				 static_cast<int>(dense->rows()),
				 static_cast<int>(dense->cols()) };
		}
		const MatrixType& d = mat.dense();
		return { const_cast<ComplexOrRealType*>(&d(0, 0)),
			 static_cast<int>(d.rows()),
			 static_cast<int>(d.cols()) };
	}

	// Needs to be public so we can use lambdas with Cuda.
public:

	struct PackTask {
		const ComplexOrRealType* src;
		ComplexOrRealType*       dst;
		int                      m, n, lds, ldd;
	};

	void setup_()
	{
		Kokkos::Profiling::ScopedRegion region("BatchedGemmKokkos::setup");

		const SizeType npatches  = initKron_.numberOfPatches(InitKronType::OLD);
		const SizeType noperator = initKron_.connections();

		// Per-patch sizes and padded leading dimensions.
		VectorSizeType lps(npatches, 0); // left patch size
		VectorSizeType rps(npatches, 0); // right patch size
		VectorSizeType ldA(npatches, 0); // padded lps -> Abatch leading dim
		VectorSizeType ldB(npatches, 0); // padded rps -> Bbatch / BXbatch leading dim
		VectorSizeType xyStart(npatches + 1, 0); // patch start offset in vin/vout

		{
			Kokkos::Profiling::ScopedRegion region(
			    "BatchedGemmKokkos::setup::patch_sizes");

			for (SizeType ip = 0; ip < npatches; ++ip) {
				const SizeType lg
				    = initKron_.patch(InitKronType::NEW, GenIjPatchType::LEFT)[ip];
				const SizeType rg
				    = initKron_.patch(InitKronType::NEW, GenIjPatchType::RIGHT)[ip];
				const int L1
				    = initKron_.lrs(InitKronType::NEW).left().partition(lg);
				const int L2
				    = initKron_.lrs(InitKronType::NEW).left().partition(lg + 1);
				const int R1
				    = initKron_.lrs(InitKronType::NEW).right().partition(rg);
				const int R2
				    = initKron_.lrs(InitKronType::NEW).right().partition(rg + 1);
				lps[ip] = static_cast<SizeType>(L2 - L1);
				rps[ip] = static_cast<SizeType>(R2 - R1);
				ldA[ip] = static_cast<SizeType>(
				    ialign_ * iceil(static_cast<int>(lps[ip]), ialign_));
				ldB[ip] = static_cast<SizeType>(
				    ialign_ * iceil(static_cast<int>(rps[ip]), ialign_));
			}

			xyStart[0] = 0;
			for (SizeType ip = 0; ip < npatches; ++ip)
				xyStart[ip + 1] = xyStart[ip] + lps[ip] * rps[ip];
		}
		// Column widths:
		//   AbatchCols[ip] = sum of lps[jp] over all non-zero (jp, k) connections for ip
		//   BbatchCols[ip] = sum of rps[jp] over all non-zero (jp, k) connections for ip
		VectorSizeType AbatchCols(npatches, 0);
		VectorSizeType BbatchCols(npatches, 0);

		{
			Kokkos::Profiling::ScopedRegion region(
			    "BatchedGemmKokkos::setup::patch_pointers");
			for (SizeType ip = 0; ip < npatches; ++ip) {
				for (SizeType jp = 0; jp < npatches; ++jp) {
					for (SizeType k = 0; k < noperator; ++k) {
						const MatrixDenseOrSparseType* Ap
						    = initKron_.xc(k)(ip, jp);
						const MatrixDenseOrSparseType* Bp
						    = initKron_.yc(k)(ip, jp);
						if (!Ap || !Bp)
							continue;
						if (Ap->isZero() || Bp->isZero())
							continue;
						AbatchCols[ip] += lps[jp];
						BbatchCols[ip] += rps[jp];
					}
				}
			}
		}
		Kokkos::Profiling::pushRegion("BatchedGemmKokkos::setup::intermediate");
		// Flat offsets: each patch ip occupies a contiguous slice.
		VectorSizeType AbatchOff(npatches + 1, 0);
		VectorSizeType BbatchOff(npatches + 1, 0);
		VectorSizeType BXbatchOff(npatches + 1, 0); // ldB rows x AbatchCols columns

		for (SizeType ip = 0; ip < npatches; ++ip) {
			AbatchOff[ip + 1]  = AbatchOff[ip] + ldA[ip] * AbatchCols[ip];
			BbatchOff[ip + 1]  = BbatchOff[ip] + ldB[ip] * BbatchCols[ip];
			BXbatchOff[ip + 1] = BXbatchOff[ip] + ldB[ip] * AbatchCols[ip];
		}

		const SizeType totalAbatch  = AbatchOff[npatches];
		const SizeType totalBbatch  = BbatchOff[npatches];
		const SizeType totalBXbatch = BXbatchOff[npatches];

		// Pass-2 split-k chunk count: patch ip's contraction dim
		// (AbatchCols[ip]) is split into ceil(AbatchCols[ip] /
		// kPass2ColChunk_) chunks (at least 1), so patches with a larger
		// contraction dimension automatically get proportionally more
		// (similarly-sized) chunks -- see kPass2ColChunk_ comment above.
		// All chunks of the same patch/tile atomically accumulate into the
		// same (m x n) region of d_vout_ (see matrixVector()), so unlike an
		// earlier version of this scheme, patch output size is NOT
		// multiplied by chunk count anywhere.
		VectorSizeType voutChunkCount(npatches, 0);
		for (SizeType ip = 0; ip < npatches; ++ip) {
			voutChunkCount[ip] = static_cast<SizeType>(std::max(
			    1, iceil(static_cast<int>(AbatchCols[ip]), kPass2ColChunk_)));
		}

		// Allocate device buffers early and create host mirrors for efficient H2D
		ExecutionSpace exec_for_alloc;
		d_flatAbatch_ = ScalarView(
		    Kokkos::view_alloc(exec_for_alloc, Kokkos::WithoutInitializing, "d_flatAbatch"),
		    totalAbatch);
		d_flatBbatch_ = ScalarView(
		    Kokkos::view_alloc(exec_for_alloc, Kokkos::WithoutInitializing, "d_flatBbatch"),
		    totalBbatch);
		d_flatBXbatch_
		    = ScalarView(Kokkos::view_alloc(
		                     exec_for_alloc, Kokkos::WithoutInitializing, "d_flatBXbatch"),
		                 totalBXbatch);

		// Create host mirrors (likely pinned) to pack data into and then do a single
		// deep_copy
		auto h_flatAbatch = Kokkos::create_mirror_view(
		    Kokkos::view_alloc(Kokkos::WithoutInitializing), d_flatAbatch_);
		auto h_flatBbatch = Kokkos::create_mirror_view(
		    Kokkos::view_alloc(Kokkos::WithoutInitializing), d_flatBbatch_);

		// Build GEMM arg lists while packing matrices.
		std::vector<GemmArgs> pass1_args;
		pass1_args.reserve(npatches * npatches * noperator);
		std::vector<GemmArgs> pass2_args;
		pass2_args.reserve(npatches);

		// Instead of performing many serial lacpy calls, collect pack tasks and
		// perform them in parallel on the host mirror.
		std::vector<PackTask> packTasks;
		packTasks.reserve(npatches * npatches); // heuristic

		Kokkos::Profiling::popRegion();

		{
			Kokkos::Profiling::ScopedRegion region("BatchedGemmKokkos::setup::gemm");

			for (SizeType ip = 0; ip < npatches; ++ip) {
				long long colA = 0; // column cursor in Abatch[ip]  (= BXbatch[ip])
				long long colB = 0; // column cursor in Bbatch[ip]

				for (SizeType jp = 0; jp < npatches; ++jp) {
					for (SizeType k = 0; k < noperator; ++k) {
						const MatrixDenseOrSparseType* Amat
						    = initKron_.xc(k)(ip, jp);
						const MatrixDenseOrSparseType* Bmat
						    = initKron_.yc(k)(ip, jp);
						if (!Amat || !Bmat)
							continue;
						if (Amat->isZero() || Bmat->isZero())
							continue;

						MatRef Aref
						    = getMatRef(*Amat); // left operator (xc)
						MatRef Bref
						    = getMatRef(*Bmat); // right operator (yc)
						assert(Aref.ptr && Bref.ptr);

						const int mA  = static_cast<int>(lps[ip]);
						const int nAk = static_cast<int>(
						    lps[jp]); // A cols = BX cols per op
						const int mB = static_cast<int>(rps[ip]);
						const int nBk
						    = static_cast<int>(rps[jp]); // B cols = X rows

						// Record pack task for A (left operator)
						PackTask tA;
						tA.src = Aref.ptr;
						tA.dst = reinterpret_cast<ComplexOrRealType*>(
						             h_flatAbatch.data())
						    + AbatchOff[ip]
						    + colA * static_cast<long long>(ldA[ip]);
						tA.m   = mA;
						tA.n   = nAk;
						tA.lds = Aref.rows;
						tA.ldd = static_cast<int>(ldA[ip]);
						packTasks.push_back(tA);

						// Record pack task for B (right operator)
						PackTask tB;
						tB.src = Bref.ptr;
						tB.dst = reinterpret_cast<ComplexOrRealType*>(
						             h_flatBbatch.data())
						    + BbatchOff[ip]
						    + colB * static_cast<long long>(ldB[ip]);
						tB.m   = mB;
						tB.n   = nBk;
						tB.lds = Bref.rows;
						tB.ldd = static_cast<int>(ldB[ip]);
						packTasks.push_back(tB);

						// Pass 1 GEMM: C(mB, nAk) = Bbatch_block(mB, nBk) *
						// X_jp(nBk, nAk)
						GemmArgs a1;
						a1.m   = mB;
						a1.n   = nAk;
						a1.k   = nBk;
						a1.lda = static_cast<int>(
						    ldB[ip]); // Bbatch leading dim
						a1.ldb = nBk; // X[jp] no padding: ld = rps[jp]
						a1.ldc = static_cast<int>(
						    ldB[ip]); // BXbatch leading dim
						a1.a_off = static_cast<long long>(BbatchOff[ip])
						    + colB * static_cast<long long>(ldB[ip]);
						a1.b_off = static_cast<long long>(xyStart[jp]);
						a1.c_off = static_cast<long long>(BXbatchOff[ip])
						    + colA * static_cast<long long>(ldB[ip]);

						// Split the m (=mB) dimension into independent
						// row-blocks: rows of C/A are not a contraction
						// dimension, so each row-block is a fully
						// self-contained GEMM and no reduction step is
						// needed (unlike the k-split used for Pass 2).
						for (int r0 = 0; r0 < a1.m; r0 += kPass1RowChunk_) {
							const int mc
							    = std::min(kPass1RowChunk_, a1.m - r0);
							GemmArgs ar = a1;
							ar.m        = mc;
							ar.a_off    = a1.a_off + r0;
							ar.c_off    = a1.c_off + r0;
							pass1_args.push_back(ar);
						}

						colA += nAk;
						colB += nBk;
					}
				}

				// Pass 2 GEMM: Y[ip](mB, mA) = BXbatch[ip](mB, k) * Abatch[ip](mA,
				// k)^T
				const int totalCols = static_cast<int>(colA); // = AbatchCols[ip]
				GemmArgs  a2;
				a2.m   = static_cast<int>(rps[ip]);
				a2.n   = static_cast<int>(lps[ip]);
				a2.k   = totalCols;
				a2.lda = static_cast<int>(ldB[ip]); // BXbatch leading dim
				a2.ldb = static_cast<int>(ldA[ip]); // Abatch leading dim
				a2.ldc
				    = static_cast<int>(rps[ip]); // Y[ip] no padding: ld = rps[ip]

				// Split the k (=totalCols) dimension into
				// voutChunkCount[ip] slices (proportional to this patch's
				// own k -- see kPass2ColChunk_ comment above), AND tile the
				// (m x n) output into at most kPass2TileDim_ x
				// kPass2TileDim_ blocks (see kPass2TileDim_ comment above,
				// bounds the Pass 2 per-team scratch buffer size). All
				// k-slices of the same (patch, output tile) write to the
				// SAME (mTile x nTile) region of d_vout_ via atomic add in
				// matrixVector(); output tiles are independent (no
				// reduction needed across tiles). Slices with 0 k-width
				// still get emitted so that TeamVectorGemmInternal's
				// beta==0 zero-fill produces a well-defined (zero) partial
				// result to add.
				const SizeType nChunks = voutChunkCount[ip];
				const int chunkCols
				    = iceil(totalCols, static_cast<int>(nChunks)); // columns per chunk
				for (int r0 = 0; r0 < a2.m; r0 += kPass2TileDim_) {
					const int mTile = std::min(kPass2TileDim_, a2.m - r0);
					for (int c0 = 0; c0 < a2.n; c0 += kPass2TileDim_) {
						const int nTile = std::min(kPass2TileDim_, a2.n - c0);

						long long colCursor = 0;
						for (SizeType c = 0; c < nChunks; ++c) {
							const int kc = std::max(0,
							    std::min(chunkCols,
							        totalCols - static_cast<int>(colCursor)));
							// Clamp the column offset so pointer arithmetic
							// never strays past this patch's allocated
							// column range (kc==0 chunks still need a valid,
							// in-bounds offset even though no data will be
							// read/written there).
							const long long colOff = std::min(
							    colCursor, static_cast<long long>(totalCols));

							GemmArgs ac = a2;
							ac.m        = mTile;
							ac.n        = nTile;
							ac.k        = kc;
							ac.a_off    = static_cast<long long>(BXbatchOff[ip])
							    + colOff * static_cast<long long>(ldB[ip]) + r0;
							ac.b_off = static_cast<long long>(AbatchOff[ip])
							    + colOff * static_cast<long long>(ldA[ip]) + c0;
							// ac.ldc keeps the FULL patch row count
							// (a2.ldc = rps[ip], unchanged from above) --
							// it is the real addressing stride into
							// d_vout_, distinct from ac.m (this tile's row
							// count) -- see matrixVector().
							ac.c_off = static_cast<long long>(xyStart[ip]) + r0
							    + static_cast<long long>(c0) * a2.ldc;
							pass2_args.push_back(ac);

							colCursor += chunkCols;
						}
					}
				}
			}

			// Execute pack tasks in parallel on the host mirror
			if (!packTasks.empty()) {
				using HostExec = Kokkos::DefaultHostExecutionSpace;

				auto      tasksPtr = packTasks.data();
				const int ntasks   = static_cast<int>(packTasks.size());
				Kokkos::parallel_for(
				    "BatchedGemmKokkos::host_pack",
				    Kokkos::RangePolicy<HostExec>(0, ntasks),
				    KOKKOS_LAMBDA(const int idx) {
					    const PackTask& t = tasksPtr[idx];
					    // column-major copy: for each column j copy m entries
					    for (int j = 0; j < t.n; ++j) {
						    std::memcpy(
						        t.dst + static_cast<long long>(j) * t.ldd,
						        t.src + static_cast<long long>(j) * t.lds,
						        static_cast<size_t>(t.m)
						            * sizeof(ComplexOrRealType));
					    }
				    });
				Kokkos::fence();
			}
		}
		{
			Kokkos::Profiling::ScopedRegion region("BatchedGemmKokkos::setup::rest");

			ExecutionSpace exec;

			// d_flatAbatch_, d_flatBbatch_, d_flatBXbatch_ were allocated above as
			// device views.
			d_vin_ = ScalarView(
			    Kokkos::view_alloc(exec, Kokkos::WithoutInitializing, "d_vin"),
			    xyStart[npatches]);
			d_vout_ = ScalarView(
			    Kokkos::view_alloc(exec, Kokkos::WithoutInitializing, "d_vout"),
			    xyStart[npatches]);

			nbatch1_ = pass1_args.size();
			nbatch2_ = pass2_args.size(); // == sum(voutChunkCount[ip])

			// Copy packed host mirrors to device in a single shot
			Kokkos::deep_copy(exec, d_flatAbatch_, h_flatAbatch);
			Kokkos::deep_copy(exec, d_flatBbatch_, h_flatBbatch);

			d_pass1_ = DevArgsView("d_pass1", nbatch1_);
			d_pass2_ = DevArgsView("d_pass2", nbatch2_);
			{
				Kokkos::View<GemmArgs*, Kokkos::HostSpace, Kokkos::MemoryUnmanaged>
				    h1(pass1_args.data(), nbatch1_);
				Kokkos::View<GemmArgs*, Kokkos::HostSpace, Kokkos::MemoryUnmanaged>
				    h2(pass2_args.data(), nbatch2_);
				Kokkos::deep_copy(exec, d_pass1_, h1);
				Kokkos::deep_copy(exec, d_pass2_, h2);
			}
			exec.fence();

			{
				PsimagLite::OstringStream msg(std::cout.precision());
				msg() << "setup done: npatches=" << npatches
				      << " noperator=" << noperator << " pass1_batches=" << nbatch1_
				      << " pass2_batches=" << nbatch2_ << " Abatch=" << totalAbatch
				      << "elems"
				      << " Bbatch=" << totalBbatch << "elems"
				      << " BXbatch=" << totalBXbatch << "elems";
				progress_.printline(msg, std::cout);
			}
		}
	}

private:

	// -----------------------------------------------------------------------

	const InitKronType&           initKron_;
	PsimagLite::ProgressIndicator progress_;
	mutable VectorMatrixType      garbage_; // owns sparse->dense expansions (freed in dtor)

	SizeType nbatch1_ = 0; // number of pass-1 GEMMs
	SizeType nbatch2_ = 0; // number of pass-2 GEMMs (== sum over patches/tiles/chunks)

	ScalarView         d_flatAbatch_; // left-operator matrices, device, persistent
	ScalarView         d_flatBbatch_; // right-operator matrices, device, persistent
	mutable ScalarView d_flatBXbatch_; // intermediate BX work buffer, per-call
	mutable ScalarView d_vin_; // input vector on device, per-call
	mutable ScalarView d_vout_; // output vector on device, per-call
	DevArgsView        d_pass1_; // pass-1 GEMM parameters, persistent
	DevArgsView        d_pass2_; // pass-2 GEMM parameters, persistent
};

} // namespace Dmrg
#endif // BATCHED_GEMM_KOKKOS_H
