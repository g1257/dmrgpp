#include "util.h"
#include <Kokkos_Core.hpp>
#include <Kokkos_Profiling_ScopedRegion.hpp>
#include <PsimagLite/KokkosType.h>

template <typename ComplexOrRealType>
void csr_to_den(const PsimagLite::CrsMatrix<ComplexOrRealType>& a,
                PsimagLite::Matrix<ComplexOrRealType>&          a_)
{
	int ia = 0;
	int ja = 0;

	int nrow_A = a.row();
	int ncol_A = a.col();

	for (ja = 0; ja < ncol_A; ja++) {
		for (ia = 0; ia < nrow_A; ia++) {
			a_(ia, ja) = 0;
		};
	};

	for (ia = 0; ia < nrow_A; ia++) {
		int istarta = a.getRowPtr(ia);
		int ienda   = a.getRowPtr(ia + 1);
		int ka      = 0;
		for (ka = istarta; ka < ienda; ka++) {
			ComplexOrRealType aij = a.getValue(ka);
			int               ja  = a.getCol(ka);
			a_(ia, ja)            = aij;
		};
	};
}

template <typename ComplexOrRealType>
void csr_kron_mult_method(const int  imethod,
                          const char transA,
                          const char transB,

                          const PsimagLite::CrsMatrix<ComplexOrRealType>& a,

                          const PsimagLite::CrsMatrix<ComplexOrRealType>& b,

                          const PsimagLite::MatrixNonOwned<const ComplexOrRealType>& yin,
                          PsimagLite::MatrixNonOwned<ComplexOrRealType>&             xout)
{
	Kokkos::Profiling::ScopedRegion region("PsimagLite::csr_kron_mult_method");

	const bool is_complex   = PsimagLite::IsComplexNumber<ComplexOrRealType>::True;
	const int  isTransA     = (transA == 'T') || (transA == 't');
	const int  isTransB     = (transB == 'T') || (transB == 't');
	const int  isConjTransA = (transA == 'C') || (transA == 'c');
	const int  isConjTransB = (transB == 'C') || (transB == 'c');

	const int nrow_A = a.rows();
	const int ncol_A = a.cols();
	const int nrow_B = b.rows();
	const int ncol_B = b.cols();

	const int nrow_1 = (isTransA || isConjTransA) ? ncol_A : nrow_A;
	const int ncol_1 = (isTransA || isConjTransA) ? nrow_A : ncol_A;
	const int nrow_2 = (isTransB || isConjTransB) ? ncol_B : nrow_B;
	const int ncol_2 = (isTransB || isConjTransB) ? nrow_B : ncol_B;

	const int nrow_X = nrow_2;
	const int ncol_X = nrow_1;
	const int nrow_Y = ncol_2;
	const int ncol_Y = ncol_1;

	assert((imethod == 1) || (imethod == 2) || (imethod == 3));

	bool no_work = (csr_is_zeros(a) || csr_is_zeros(b));
	if (no_work) {
		return;
	};
	/*
	 *   -------------------------------------------------------------
	 *   A and B in compressed sparse ROW format
	 *
	 *   X += kron( op(A), op(B)) * Y
	 *   X +=  op(B) * Y * transpose(op(A))
	 *
	 *   nrow_X = nrow_2,   ncol_X = nrow_1
	 *   nrow_Y = ncol_2,   nrow_Y = nrow_2
	 *
	 *   that can be computed as either
	 *   imethod == 1
	 *
	 *   X(ib,ia) +=  (B(ib,jb) * Y( jb,ja)) * transpose( A(ia,ja))
	 *
	 *   X(ix,jx) +=  (B2(ib2,jb2) * Y(iy,jy)) * transpose(A1(ia1,ja1))
	 *
	 *
	 *
	 *   X(ib,ia) += (B(ib,jb) * Y(jb,ja) ) * transpose(A(ia,ja)   or
	 *               BY(ib,ja) = B(ib,jb)*Y(jb,ja)
	 *               BY is nrow_B by ncol_A, need   2*nnz(B)*ncolA flops
	 *
	 *   X(ib,ia) +=   BY(ib,ja) * transpose(A(ia,ja)) need 2*nnz(A)*nrowB flops
	 *
	 *   imethod == 2
	 *
	 *   X(ib,ia) += B(ib,jb) * (Y(jb,ja) * transpose(A))    or
	 *                YAt(jb,ia) = Y(jb,ja) * transpose(A(ia,ja))
	 *                YAt is ncolB by nrowA, need 2*nnz(A) * ncolB flops
	 *
	 *   X(ib,ia) += B(ib,jb) * YAt(jb,ia)  need nnz(B) * nrowA flops
	 *
	 *   imethod == 3
	 *
	 *   X += kron(A,B) * Y   by visiting all non-zero entries in A, B
	 *
	 *   this is feasible only if A and B are very sparse, need nnz(A)*nnz(B) flops
	 *   -------------------------------------------------------------
	 */

	if (imethod == 1) {
		Kokkos::Profiling::ScopedRegion region("PsimgLite::csr_kron_mult_method::imethod1");

		/*
		 *  --------------------------------------------
		 *  BY(ib,ja) = (B(ib,jb))*Y(jb,ja)
		 *
		 *  X(ib,ia) += BY(ib,ja ) * transpose(A(ia,ja))
		 *  --------------------------------------------
		 */

		int                                                 nrow_BY = nrow_X;
		int                                                 ncol_BY = ncol_Y;
		PsimagLite::Matrix<ComplexOrRealType>               by_(nrow_BY, ncol_BY);
		PsimagLite::MatrixNonOwned<ComplexOrRealType>       byRef(by_);
		PsimagLite::MatrixNonOwned<const ComplexOrRealType> byConstRef(by_);
		/*
		 * ---------------
		 * setup BY
		 * ---------------
		 */

		{
			int iby = 0;
			int jby = 0;

			// not needed FIXME
			for (jby = 0; jby < ncol_BY; jby++) {
				for (iby = 0; iby < nrow_BY; iby++) {
					by_(iby, jby) = 0;
				};
			};
		}

		{
			/*
			 * ------------------------------
			 * BY(ib,ja)  = B(ib,jb)*Y(jb,ja)
			 * ------------------------------
			 */
			// const char trans = (isTransB) ? 'T' : 'N';
			const char trans = transB;
			csr_matmul_pre(trans,
			               b,

			               nrow_Y,
			               ncol_Y,
			               yin,

			               nrow_BY,
			               ncol_BY,
			               byRef);
		}

		{
			/*
			 * -------------------------------------------
			 * X(ib,ia) += BY(ib,ja) * transpose(A(ia,ja))
			 * -------------------------------------------
			 */
			/*
			 * ---------------------------------
			 * note trans = 'Z' mean use conj(A)
			 * ---------------------------------
			 */
			const char trans = isTransA ? 'N' : (isConjTransA ? 'Z' : 'T');
			csr_matmul_post(trans,
			                a,

			                nrow_BY,
			                ncol_BY,
			                byConstRef,

			                nrow_X,
			                ncol_X,
			                xout);
		}
	} else if (imethod == 2) {
		Kokkos::Profiling::ScopedRegion region("PsimgLite::csr_kron_mult_method::imethod2");

		/*
		 * ---------------------
		 * YAt(jb,ia) = Y(jb,ja) * tranpose(A(ia,ja))
		 * X(ib,ia) += B(ib,jb) * YAt(jb,ia)
		 * ---------------------
		 */

		int                                                 nrow_YAt = nrow_Y;
		int                                                 ncol_YAt = ncol_X;
		PsimagLite::Matrix<ComplexOrRealType>               yat_(nrow_YAt, ncol_YAt);
		PsimagLite::MatrixNonOwned<ComplexOrRealType>       yatRef(yat_);
		PsimagLite::MatrixNonOwned<const ComplexOrRealType> yatConstRef(yat_);

		/*
		 * ----------------
		 * setup YAt(jb,ia)
		 * ----------------
		 */

		{
			int iy = 0;
			int jy = 0;

			// not needed FIXME
			for (jy = 0; jy < ncol_YAt; jy++) {
				for (iy = 0; iy < nrow_YAt; iy++) {
					yat_(iy, jy) = 0;
				};
			};
		}

		{
			/*
			 * ---------------------
			 * YAt(jb,ia) = Y(jb,ja) * tranpose(A(ia,ja)
			 * ---------------------
			 */
			/*
			 * ---------------------------------
			 * note trans = 'Z' mean use conj(A)
			 * ---------------------------------
			 */
			const char transa = isTransA ? 'N' : (isConjTransA ? 'Z' : 'T');
			csr_matmul_post(transa,
			                a,

			                nrow_Y,
			                ncol_Y,
			                yin,

			                nrow_YAt,
			                ncol_YAt,
			                yatRef);
		}

		{
			/*
			 * ------------
			 * X(ib,ia) += B(ib,jb) * YAt(jb,ia)
			 * ------------
			 */

			// const char trans = (isTransB) ? 'T' : 'N';
			const char trans = transB;

			csr_matmul_pre(trans,
			               b,

			               nrow_YAt,
			               ncol_YAt,
			               yatConstRef,

			               nrow_X,
			               ncol_X,
			               xout);
		}
	} else if (imethod == 3) {
		Kokkos::Profiling::ScopedRegion region("PsimgLite::csr_kron_mult_method::imethod3");
		/*
		 * ---------------------------------------------
		 * C = kron(A,B)
		 * C([ib,ia], [jb,ja]) = A(ia,ja)*B(ib,jb)
		 * X([ib,ia]) += C([ib,ia],[jb,ja]) * Y([jb,ja])
		 * ---------------------------------------------
		 */

		using ExecutionSpace = Kokkos::DefaultExecutionSpace;
		using KokkosScalar   = typename PsimagLite::KokkosType<ComplexOrRealType>::type;

		const int nnzA = a.nonZeros();
		const int nnzB = b.nonZeros();

		/*
		 * CPU fast-path: for small problems, direct loops avoid GPU transfer
		 * overhead (~400+ µs/call on discrete GPUs). The threshold below (tunable)
		 * targets the crossover where GPU kernel savings outweigh copy latency.
		 */
		const size_t totalPairs  = static_cast<size_t>(nnzA) * static_cast<size_t>(nnzB);
		static constexpr size_t kGpuThreshold = 100000;
		if (totalPairs < kGpuThreshold) {
			Kokkos::Profiling::ScopedRegion cpuRegion("PsimgLite::csr_kron_mult_method::imethod3::cpu");
			for (int ia = 0; ia < nrow_A; ++ia) {
				const int istart_a = a.getRowPtr(ia);
				const int iend_a   = a.getRowPtr(ia + 1);
				for (int ka = istart_a; ka < iend_a; ++ka) {
					const int         ja  = a.getCol(ka);
					ComplexOrRealType aij = a.getValue(ka);
					if constexpr (is_complex)
						if (isConjTransA)
							aij = PsimagLite::conj(aij);
					for (int ib = 0; ib < nrow_B; ++ib) {
						const int istart_b = b.getRowPtr(ib);
						const int iend_b   = b.getRowPtr(ib + 1);
						for (int kb = istart_b; kb < iend_b; ++kb) {
							const int         jb  = b.getCol(kb);
							ComplexOrRealType bij = b.getValue(kb);
							if constexpr (is_complex)
								if (isConjTransB)
									bij = PsimagLite::conj(bij);

							const int ix = (isTransB || isConjTransB) ? jb : ib;
							const int jx = (isTransA || isConjTransA) ? ja : ia;
							const int iy = (isTransB || isConjTransB) ? ib : jb;
							const int jy = (isTransA || isConjTransA) ? ia : ja;

							xout(ix, jx) += aij * bij * yin(iy, jy);
						}
					}
				}
			}
			return;
		}

		/*
		 * GPU path for large problems.
		 *
		 * Use compact rowptr arrays instead of expanded row-index arrays.
		 * All four are unmanaged host views wrapping the existing CrsMatrix
		 * storage (zero extra allocation on the host side).
		 *
		 * Kernel structure:  TeamPolicy(nrow_A, AUTO)
		 *   league rank  → A row ia  (one team per A row)
		 *   team threads → B row ib  (one thread per B row)
		 *   serial inner → A / B nonzeros within (ia, ib)
		 *
		 * For the non-transpose case this layout eliminates atomic_add:
		 *   jx = ia   — unique per team   → no inter-team write conflict
		 *   ix = ib   — unique per thread → no intra-team write conflict
		 * x_dev(ib, ia) is written by exactly one (team, thread) pair.
		 *
		 * With LayoutLeft both the x_dev writes and the y_dev reads are
		 * stride-1 (ib / jb vary across threads, ia / ja fixed per team).
		 */
		Kokkos::View<const int*, Kokkos::HostSpace, Kokkos::MemoryUnmanaged>
		    A_rowptr_h(&a.getRowPtr(0), nrow_A + 1);
		Kokkos::View<const int*, Kokkos::HostSpace, Kokkos::MemoryUnmanaged>
		    A_col_h(&a.getCol(0), nnzA);
		Kokkos::View<const KokkosScalar*, Kokkos::HostSpace, Kokkos::MemoryUnmanaged>
		    A_val_h(reinterpret_cast<const KokkosScalar*>(&a.getValue(0)), nnzA);

		Kokkos::View<const int*, Kokkos::HostSpace, Kokkos::MemoryUnmanaged>
		    B_rowptr_h(&b.getRowPtr(0), nrow_B + 1);
		Kokkos::View<const int*, Kokkos::HostSpace, Kokkos::MemoryUnmanaged>
		    B_col_h(&b.getCol(0), nnzB);
		Kokkos::View<const KokkosScalar*, Kokkos::HostSpace, Kokkos::MemoryUnmanaged>
		    B_val_h(reinterpret_cast<const KokkosScalar*>(&b.getValue(0)), nnzB);

		auto A_rowptr_dev = Kokkos::create_mirror_view_and_copy(ExecutionSpace{}, A_rowptr_h);
		auto A_col_dev    = Kokkos::create_mirror_view_and_copy(ExecutionSpace{}, A_col_h);
		auto A_val_dev    = Kokkos::create_mirror_view_and_copy(ExecutionSpace{}, A_val_h);
		auto B_rowptr_dev = Kokkos::create_mirror_view_and_copy(ExecutionSpace{}, B_rowptr_h);
		auto B_col_dev    = Kokkos::create_mirror_view_and_copy(ExecutionSpace{}, B_col_h);
		auto B_val_dev    = Kokkos::create_mirror_view_and_copy(ExecutionSpace{}, B_val_h);

		auto yin_host = Kokkos::View<const KokkosScalar**,
		                             Kokkos::LayoutLeft,
		                             Kokkos::HostSpace,
		                             Kokkos::MemoryUnmanaged>(
		    reinterpret_cast<const KokkosScalar*>(&yin(0, 0)), nrow_Y, ncol_Y);
		auto y_dev = Kokkos::create_mirror_view_and_copy(ExecutionSpace{}, yin_host);

		auto x_dev = Kokkos::View<KokkosScalar**, Kokkos::LayoutLeft>("x_dev", nrow_X, ncol_X);
		Kokkos::deep_copy(x_dev, KokkosScalar(0));

		using TeamPolicy = Kokkos::TeamPolicy<ExecutionSpace>;
		using TeamMember = typename TeamPolicy::member_type;

		if (!isTransA && !isConjTransA && !isTransB && !isConjTransB) {
			/*
			 * Non-transpose path: ix = ib (unique per thread), jx = ia (unique
			 * per team).  Each (ia, ib) cell of x_dev is touched by exactly one
			 * (team, thread) pair → plain += without atomic_add is correct.
			 */
			Kokkos::parallel_for(
			    "csr_kron_mult::imethod3_nn",
			    TeamPolicy(nrow_A, Kokkos::AUTO),
			    KOKKOS_LAMBDA(const TeamMember& team) {
				    const int ia       = team.league_rank();
				    const int ka_begin = A_rowptr_dev(ia);
				    const int ka_end   = A_rowptr_dev(ia + 1);
				    Kokkos::parallel_for(
				        Kokkos::TeamThreadRange(team, nrow_B),
				        [=](int ib) {
					        const int    kb_begin = B_rowptr_dev(ib);
					        const int    kb_end   = B_rowptr_dev(ib + 1);
					        KokkosScalar acc      = 0;
					        for (int ka = ka_begin; ka < ka_end; ++ka) {
						        const int    ja  = A_col_dev(ka);
						        KokkosScalar aij = A_val_dev(ka);
						        for (int kb = kb_begin; kb < kb_end; ++kb)
							        acc += aij * B_val_dev(kb) * y_dev(B_col_dev(kb), ja);
					        }
					        x_dev(ib, ia) += acc;
				        });
			    });
		} else {
			/*
			 * General (transpose / conjugate) path: ix or jx may not be unique
			 * per (team, thread) pair → atomic_add required.
			 */
			Kokkos::parallel_for(
			    "csr_kron_mult::imethod3_gen",
			    TeamPolicy(nrow_A, Kokkos::AUTO),
			    KOKKOS_LAMBDA(const TeamMember& team) {
				    const int ia       = team.league_rank();
				    const int ka_begin = A_rowptr_dev(ia);
				    const int ka_end   = A_rowptr_dev(ia + 1);
				    Kokkos::parallel_for(
				        Kokkos::TeamThreadRange(team, nrow_B),
				        [=](int ib) {
					        const int kb_begin = B_rowptr_dev(ib);
					        const int kb_end   = B_rowptr_dev(ib + 1);
					        for (int ka = ka_begin; ka < ka_end; ++ka) {
						        const int    ja  = A_col_dev(ka);
						        KokkosScalar aij = A_val_dev(ka);
						        if constexpr (is_complex)
							        if (isConjTransA)
								        aij = Kokkos::conj(aij);
						        for (int kb = kb_begin; kb < kb_end; ++kb) {
							        const int    jb  = B_col_dev(kb);
							        KokkosScalar bij = B_val_dev(kb);
							        if constexpr (is_complex)
								        if (isConjTransB)
									        bij = Kokkos::conj(bij);
							        const int ix = (isTransB || isConjTransB) ? jb : ib;
							        const int jx = (isTransA || isConjTransA) ? ja : ia;
							        const int iy = (isTransB || isConjTransB) ? ib : jb;
							        const int jy = (isTransA || isConjTransA) ? ia : ja;
							        Kokkos::atomic_add(&x_dev(ix, jx), aij * bij * y_dev(iy, jy));
						        }
					        }
				        });
			    });
		}

		// copy result back and accumulate into xout
		auto xhost = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace{}, x_dev);
		for (int ix = 0; ix < nrow_X; ++ix)
			for (int jx = 0; jx < ncol_X; ++jx)
				xout(ix, jx) += static_cast<ComplexOrRealType>(xhost(ix, jx));
	};
}

template <typename ComplexOrRealType>
void csr_kron_mult_method(const int                                                   imethod,
                          const char                                                  transA,
                          const char                                                  transB,
                          const PsimagLite::CrsMatrix<ComplexOrRealType>&             a,
                          const PsimagLite::CrsMatrix<ComplexOrRealType>&             b,
                          const typename PsimagLite::Vector<ComplexOrRealType>::Type& yin_,
                          SizeType                                                    offsetY,
                          typename PsimagLite::Vector<ComplexOrRealType>::Type&       xout_,
                          SizeType                                                    offsetX)

{
	const int isTransA     = (transA == 'T') || (transA == 't');
	const int isTransB     = (transB == 'T') || (transB == 't');
	const int isConjTransA = (transA == 'C') || (transA == 'c');
	const int isConjTransB = (transB == 'C') || (transB == 'c');

	const int nrow_A = a.rows();
	const int ncol_A = a.cols();
	const int nrow_B = b.rows();
	const int ncol_B = b.cols();

	const int nrow_1 = (isTransA || isConjTransA) ? ncol_A : nrow_A;
	const int ncol_1 = (isTransA || isConjTransA) ? nrow_A : ncol_A;
	const int nrow_2 = (isTransB || isConjTransB) ? ncol_B : nrow_B;
	const int ncol_2 = (isTransB || isConjTransB) ? nrow_B : ncol_B;

	const int                                           nrow_X = nrow_2;
	const int                                           ncol_X = nrow_1;
	const int                                           nrow_Y = ncol_2;
	const int                                           ncol_Y = ncol_1;
	PsimagLite::MatrixNonOwned<const ComplexOrRealType> yin(nrow_Y, ncol_Y, yin_, offsetY);
	PsimagLite::MatrixNonOwned<ComplexOrRealType>       xout(nrow_X, ncol_X, xout_, offsetX);
	csr_kron_mult_method(imethod, transA, transB, a, b, yin, xout);
}

template <typename ComplexOrRealType>
void csr_kron_mult(const char                                                  transA,
                   const char                                                  transB,
                   const PsimagLite::CrsMatrix<ComplexOrRealType>&             a,
                   const PsimagLite::CrsMatrix<ComplexOrRealType>&             b,
                   const typename PsimagLite::Vector<ComplexOrRealType>::Type& yin,
                   SizeType                                                    offsetY,
                   typename PsimagLite::Vector<ComplexOrRealType>::Type&       xout,
                   SizeType                                                    offsetX,
                   const typename PsimagLite::Real<ComplexOrRealType>::Type    denseFlopDiscount)
{
	/*
	 *   -------------------------------------------------------------
	 *   A and B in compressed sparse ROW format
	 *
	 *   X += kron( A, B) * Y
	 *   that can be computed as either
	 *   imethod == 1
	 *
	 *   X(ib,ia) += (B(ib,jb) * Y(jb,ja) ) * transpose(A(ia,ja)   or
	 *               BY(ib,ja) = B(ib,jb)*Y(jb,ja)
	 *               BY is nrow_B by ncol_A, need   2*nnz(B)*ncolA flops
	 *
	 *   X(ib,ia) +=   BY(ib,ja) * transpose(A(ia,ja)) need 2*nnz(A)*nrowB flops
	 *
	 *   imethod == 2
	 *
	 *   X(ib,ia) += B(ib,jb) * (Y(jb,ja) * transpose(A))    or
	 *                YAt(jb,ia) = Y(jb,ja) * transpose(A(ia,ja))
	 *                YAt is ncolB by nrowA, need 2*nnz(A) * ncolB flops
	 *
	 *   X(ib,ia) += B(ib,jb) * YAt(jb,ia)  need nnz(B) * nrowA flops
	 *
	 *   imethod == 3
	 *
	 *   X += kron(A,B) * Y   by visiting all non-zero entries in A, B
	 *
	 *   this is feasible only if A and B are very sparse, need nnz(A)*nnz(B) flops
	 *   -------------------------------------------------------------
	 */
	int nnz_A = csr_nnz(a);
	int nnz_B = csr_nnz(b);

	bool no_work = (csr_is_zeros(a) || csr_is_zeros(b));
	if (no_work) {
		return;
	};

	ComplexOrRealType kron_nnz   = 0;
	ComplexOrRealType kron_flops = 0;
	int               imethod    = 1;

	const int isTransA     = (transA == 'T') || (transA == 't');
	const int isTransB     = (transB == 'T') || (transB == 't');
	const int isConjTransA = (transA == 'C') || (transA == 'c');
	const int isConjTransB = (transB == 'C') || (transB == 'c');

	const int nrow_A = a.rows();
	const int ncol_A = a.cols();
	const int nrow_B = b.rows();
	const int ncol_B = b.cols();

	// -----------------------------------
	// both A and B are considered sparse
	// -----------------------------------

	const int nrow_1 = (isTransA || isConjTransA) ? ncol_A : nrow_A;
	const int ncol_1 = (isTransA || isConjTransA) ? nrow_A : ncol_A;

	const int nrow_2 = (isTransB || isConjTransB) ? ncol_B : nrow_B;
	const int ncol_2 = (isTransB || isConjTransB) ? nrow_B : ncol_B;

	estimate_kron_cost(nrow_1,
	                   ncol_1,
	                   nnz_A,
	                   nrow_2,
	                   ncol_2,
	                   nnz_B,
	                   &kron_nnz,
	                   &kron_flops,
	                   &imethod,
	                   denseFlopDiscount);

	csr_kron_mult_method(imethod, transA, transB, a, b, yin, offsetY, xout, offsetX);
}
