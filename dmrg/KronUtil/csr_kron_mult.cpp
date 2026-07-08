#include "util.h"

#include "PsimagLiteConfig.h"
#ifdef PSIMAGLITE_USE_KOKKOS
#include <Kokkos_Core.hpp>
#include <Kokkos_Profiling_ScopedRegion.hpp>

template <typename T> struct KokkosType {
	using type = T;
};

template <typename T>
        requires(!std::is_floating_point_v<T>)
struct KokkosType<T> {
	using type = Kokkos::complex<typename T::value_type>;
};
#endif

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
void csr_kron_mult_method(const int                                                  imethod,
                          const char                                                 transA,
                          const char                                                 transB,
                          const PsimagLite::CrsMatrix<ComplexOrRealType>&            a,
                          const PsimagLite::CrsMatrix<ComplexOrRealType>&            b,
                          const PsimagLite::MatrixNonOwned<const ComplexOrRealType>& yin,
                          PsimagLite::MatrixNonOwned<ComplexOrRealType>&             xout)
{
	constexpr bool is_complex   = PsimagLite::IsComplexNumber<ComplexOrRealType>::True;
	const int      isTransA     = (transA == 'T') || (transA == 't');
	const int      isTransB     = (transB == 'T') || (transB == 't');
	const int      isConjTransA = (transA == 'C') || (transA == 'c');
	const int      isConjTransB = (transB == 'C') || (transB == 'c');

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
#ifdef PSIMAGLITE_USE_KOKKOS
		Kokkos::Profiling::ScopedRegion region("PsimgLite::csr_kron_mult_method::imethod1");
#endif
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
#ifdef PSIMAGLITE_USE_KOKKOS
		Kokkos::Profiling::ScopedRegion region("PsimgLite::csr_kron_mult_method::imethod2");
#endif
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
#ifdef PSIMAGLITE_USE_KOKKOS
		Kokkos::Profiling::ScopedRegion region("PsimgLite::csr_kron_mult_method::imethod3");
#endif
		/*
		 * ---------------------------------------------
		 * C = kron(A,B)
		 * C([ib,ia], [jb,ja]) = A(ia,ja)*B(ib,jb)
		 * X([ib,ia]) += C([ib,ia],[jb,ja]) * Y([jb,ja])
		 * ---------------------------------------------
		 */

#ifdef PSIMAGLITE_USE_KOKKOS
		bool doTransA = (isTransA || isConjTransA);
		bool doTransB = (isTransB || isConjTransB);

		// Build arrays for A
		int nnzA = 0;
		for (int ia1 = 0; ia1 < nrow_A; ++ia1)
			nnzA += (a.getRowPtr(ia1 + 1) - a.getRowPtr(ia1));
		std::vector<int> a_rowptr(nrow_A + 1);
		for (int i = 0; i <= nrow_A; ++i)
			a_rowptr[i] = a.getRowPtr(i);
		std::vector<int>               a_cols(nnzA);
		std::vector<ComplexOrRealType> a_vals(nnzA);
		for (int k = 0, ia1 = 0; ia1 < nrow_A; ++ia1) {
			int ist  = a_rowptr[ia1];
			int iend = a_rowptr[ia1 + 1];
			for (int kk = ist; kk < iend; ++kk, ++k) {
				a_cols[k] = a.getCol(kk);
				a_vals[k] = a.getValue(kk);
			}
		}

		// Build arrays for B
		int nnzB = 0;
		for (int ib1 = 0; ib1 < nrow_B; ++ib1)
			nnzB += (b.getRowPtr(ib1 + 1) - b.getRowPtr(ib1));
		std::vector<int> b_rowptr(nrow_B + 1);
		for (int i = 0; i <= nrow_B; ++i)
			b_rowptr[i] = b.getRowPtr(i);
		std::vector<int>               b_cols(nnzB);
		std::vector<ComplexOrRealType> b_vals(nnzB);
		std::vector<int>               b_row_of(nnzB);
		for (int k = 0, ib1 = 0; ib1 < nrow_B; ++ib1) {
			int ist  = b_rowptr[ib1];
			int iend = b_rowptr[ib1 + 1];
			for (int kk = ist; kk < iend; ++kk, ++k) {
				b_cols[k]   = b.getCol(kk);
				b_vals[k]   = b.getValue(kk);
				b_row_of[k] = ib1; // remember the row for this nonzero
			}
		}

		// Build CSC-like index for B: group nonzeros by column jb
		std::vector<int> colCountsB(ncol_B, 0);
		for (int k = 0; k < nnzB; ++k) {
			int jb = b_cols[k];
			if (jb >= 0 && jb < ncol_B)
				++colCountsB[jb];
		}
		std::vector<int> colPtrB(ncol_B + 1, 0);
		for (int i = 0; i < ncol_B; ++i)
			colPtrB[i + 1] = colPtrB[i] + colCountsB[i];
		std::vector<int> idxListB(nnzB);
		{
			std::vector<int> curB(colPtrB.begin(), colPtrB.end());
			for (int k = 0; k < nnzB; ++k) {
				int jb        = b_cols[k];
				int pos       = curB[jb]++;
				idxListB[pos] = k; // store index into b_* arrays
			}
		}

		// Map A nonzeros to target columns (jx)
		int              nTarget = ncol_X;
		std::vector<int> colCounts(nTarget, 0);
		std::vector<int> rindexA(nnzA);
		for (int ia1 = 0, k = 0; ia1 < nrow_A; ++ia1) {
			int ist  = a_rowptr[ia1];
			int iend = a_rowptr[ia1 + 1];
			for (int kk = ist; kk < iend; ++kk, ++k) {
				rindexA[k] = ia1;
				int ja     = a_cols[k];
				int tgt    = doTransA ? ja : ia1;
				if (tgt >= 0 && tgt < nTarget)
					++colCounts[tgt];
			}
		}
		std::vector<int> colPtr(nTarget + 1, 0);
		for (int i = 0; i < nTarget; ++i)
			colPtr[i + 1] = colPtr[i] + colCounts[i];
		std::vector<int> idxList(nnzA);
		std::vector<int> cur(colPtr.begin(), colPtr.end());
		for (int k = 0; k < nnzA; ++k) {
			int tgt      = doTransA ? a_cols[k] : rindexA[k];
			int pos      = cur[tgt]++;
			idxList[pos] = k;
		}

		using HostExec = Kokkos::DefaultExecutionSpace;
		HostExec exec;
		Kokkos::fence();

		// Prepare output flat buffer to collect per-column results (flattened by column)
		std::vector<ComplexOrRealType> results_flat((size_t)nTarget * (size_t)nrow_X);
		// initialize from current xout
		for (int jx = 0; jx < nTarget; ++jx)
			for (int ix = 0; ix < nrow_X; ++ix)
				results_flat[(size_t)jx * (size_t)nrow_X + (size_t)ix]
				    = xout(ix, jx);
		ComplexOrRealType* results_ptr = results_flat.data();

		// prepare flat yin array to avoid capturing MatrixNonOwned
		std::vector<ComplexOrRealType> yin_flat((size_t)nrow_Y * (size_t)ncol_Y);
		for (int jy = 0; jy < ncol_Y; ++jy)
			for (int iy = 0; iy < nrow_Y; ++iy)
				yin_flat[(size_t)iy + (size_t)jy * (size_t)nrow_Y] = yin(iy, jy);
		ComplexOrRealType* yin_ptr = yin_flat.data();

		// expose raw pointers for arrays so lambda copy is cheap
		int*               a_cols_ptr   = a_cols.data();
		ComplexOrRealType* a_vals_ptr   = a_vals.data();
		int*               rindexA_ptr  = rindexA.data();
		int*               colPtr_ptr   = colPtr.data();
		int*               idxList_ptr  = idxList.data();
		int*               b_rowptr_ptr = b_rowptr.data();
		int*               b_cols_ptr   = b_cols.data();
		ComplexOrRealType* b_vals_ptr   = b_vals.data();

		// Create device Views for arrays
		using Exec     = Kokkos::DefaultExecutionSpace;
		using MemSpace = typename Exec::memory_space;

		using KokkosScalar = KokkosType<ComplexOrRealType>::type;

		Kokkos::View<int*, MemSpace> a_cols_dev(
		    Kokkos::view_alloc("a_cols_dev", Kokkos::WithoutInitializing), nnzA);
		Kokkos::View<KokkosScalar*, MemSpace> a_vals_dev(
		    Kokkos::view_alloc("a_vals_dev", Kokkos::WithoutInitializing), nnzA);
		Kokkos::View<int*, MemSpace> rindexA_dev(
		    Kokkos::view_alloc("rindexA_dev", Kokkos::WithoutInitializing), nnzA);
		Kokkos::View<int*, MemSpace> colPtr_dev(
		    Kokkos::view_alloc("colPtr_dev", Kokkos::WithoutInitializing), nTarget + 1);
		Kokkos::View<int*, MemSpace> idxList_dev(
		    Kokkos::view_alloc("idxList_dev", Kokkos::WithoutInitializing), nnzA);
		Kokkos::View<int*, MemSpace> b_rowptr_dev(
		    Kokkos::view_alloc("b_rowptr_dev", Kokkos::WithoutInitializing), nrow_B + 1);
		Kokkos::View<int*, MemSpace> b_cols_dev(
		    Kokkos::view_alloc("b_cols_dev", Kokkos::WithoutInitializing), nnzB);
		Kokkos::View<KokkosScalar*, MemSpace> b_vals_dev(
		    Kokkos::view_alloc("b_vals_dev", Kokkos::WithoutInitializing), nnzB);
		Kokkos::View<int*, MemSpace> b_row_of_dev(
		    Kokkos::view_alloc("b_row_of_dev", Kokkos::WithoutInitializing), nnzB);
		Kokkos::View<int*, MemSpace> colPtrB_dev(
		    Kokkos::view_alloc("colPtrB_dev", Kokkos::WithoutInitializing), ncol_B + 1);
		Kokkos::View<int*, MemSpace> idxListB_dev(
		    Kokkos::view_alloc("idxListB_dev", Kokkos::WithoutInitializing), nnzB);
		Kokkos::View<KokkosScalar*, MemSpace> results_dev(
		    Kokkos::view_alloc("results_dev", Kokkos::WithoutInitializing),
		    (size_t)nTarget * (size_t)nrow_X);
		Kokkos::View<KokkosScalar**, Kokkos::LayoutLeft, MemSpace> yin_dev(
		    Kokkos::view_alloc("yin_dev", Kokkos::WithoutInitializing), nrow_Y, ncol_Y);

		// fill host mirrors then deep_copy to device
		auto a_cols_h   = Kokkos::create_mirror_view(a_cols_dev);
		auto a_vals_h   = Kokkos::create_mirror_view(a_vals_dev);
		auto rindexA_h  = Kokkos::create_mirror_view(rindexA_dev);
		auto colPtr_h   = Kokkos::create_mirror_view(colPtr_dev);
		auto idxList_h  = Kokkos::create_mirror_view(idxList_dev);
		auto b_rowptr_h = Kokkos::create_mirror_view(b_rowptr_dev);
		auto b_cols_h   = Kokkos::create_mirror_view(b_cols_dev);
		auto b_vals_h   = Kokkos::create_mirror_view(b_vals_dev);
		auto b_row_of_h = Kokkos::create_mirror_view(b_row_of_dev);
		auto colPtrB_h  = Kokkos::create_mirror_view(colPtrB_dev);
		auto idxListB_h = Kokkos::create_mirror_view(idxListB_dev);
		auto results_h  = Kokkos::create_mirror_view(results_dev);
		auto yin_h      = Kokkos::create_mirror_view(yin_dev);

		for (int k = 0; k < nnzA; ++k) {
			a_cols_h(k)  = a_cols[k];
			a_vals_h(k)  = a_vals[k];
			rindexA_h(k) = rindexA[k];
			idxList_h(k) = idxList[k];
		}
		for (int i = 0; i <= nTarget; ++i)
			colPtr_h(i) = colPtr[i];
		for (int k = 0; k < nnzB; ++k) {
			b_cols_h(k)   = b_cols[k];
			b_vals_h(k)   = b_vals[k];
			b_row_of_h(k) = b_row_of[k];
		}
		for (int i = 0; i <= nrow_B; ++i)
			b_rowptr_h(i) = b_rowptr[i];
		for (int i = 0; i <= ncol_B; ++i)
			colPtrB_h(i) = colPtrB[i];
		for (int k = 0; k < nnzB; ++k)
			idxListB_h(k) = idxListB[k];
		// initialize results_h from results_ptr host data
		for (size_t i = 0; i < (size_t)nTarget * (size_t)nrow_X; ++i)
			results_h(i) = results_ptr[i];
		for (int jy = 0; jy < ncol_Y; ++jy)
			for (int iy = 0; iy < nrow_Y; ++iy)
				yin_h(iy, jy) = yin_ptr[(size_t)iy + (size_t)jy * (size_t)nrow_Y];

		Kokkos::deep_copy(a_cols_dev, a_cols_h);
		Kokkos::deep_copy(a_vals_dev, a_vals_h);
		Kokkos::deep_copy(rindexA_dev, rindexA_h);
		Kokkos::deep_copy(colPtr_dev, colPtr_h);
		Kokkos::deep_copy(idxList_dev, idxList_h);
		Kokkos::deep_copy(b_rowptr_dev, b_rowptr_h);
		Kokkos::deep_copy(b_cols_dev, b_cols_h);
		Kokkos::deep_copy(b_vals_dev, b_vals_h);
		Kokkos::deep_copy(b_row_of_dev, b_row_of_h);
		Kokkos::deep_copy(colPtrB_dev, colPtrB_h);
		Kokkos::deep_copy(idxListB_dev, idxListB_h);
		Kokkos::deep_copy(b_rowptr_dev, b_rowptr_h);
		Kokkos::deep_copy(b_cols_dev, b_cols_h);
		Kokkos::deep_copy(b_vals_dev, b_vals_h);
		Kokkos::deep_copy(results_dev, results_h);
		Kokkos::deep_copy(yin_dev, yin_h);

		// TeamPolicy over target columns; each team writes its slice of results_dev
		using team_policy = Kokkos::TeamPolicy<Exec>;
		using member_type = team_policy::member_type;
		team_policy policy(nTarget, Kokkos::AUTO);
		Kokkos::parallel_for(
		    "imethod3_kron_team", policy, KOKKOS_LAMBDA(const member_type& member) {
			    int jx = member.league_rank();
			    // Parallelize over output rows ix using team threads; selectively
			    // vectorize inner loops
			    Kokkos::parallel_for(
			        Kokkos::TeamThreadRange(member, nrow_X),
			        [&](const int ix)
			        {
				        auto acc
				            = results_dev((size_t)jx * (size_t)nrow_X + (size_t)ix);

				        int start = colPtr_dev(jx);
				        int end   = colPtr_dev(jx + 1);

				        for (int p = start; p < end; ++p) {
					        int  ka  = idxList_dev(p);
					        int  ia1 = rindexA_dev(ka);
					        int  ja  = a_cols_dev(ka);
					        auto aij = a_vals_dev(ka);
					        if constexpr (is_complex)
						        if (isConjTransA)
							        aij = Kokkos::conj(aij);

					        if (!doTransB) {
						        int ib1    = ix; // ix corresponds to ib1
						        int bstart = b_rowptr_dev(ib1);
						        int bend   = b_rowptr_dev(ib1 + 1);
						        int len    = bend - bstart;
						        if (len > 32) {
							        Kokkos::parallel_for(
							            Kokkos::ThreadVectorRange(
							                member, len),
							            [&](const int t)
							            {
								            int kb = bstart + t;
								            int jb = b_cols_dev(kb);
								            auto bij
								                = b_vals_dev(kb);
								            if constexpr (
								                is_complex)
									            if (isConjTransB)
										            bij = Kokkos::
										                conj(
										                    bij);
								            auto cij = aij * bij;
								            int  iy  = jb;
								            int  jy = doTransA ? ia1
								                               : ja;
								            acc += cij
								                * yin_dev(iy, jy);
							            });
						        } else {
							        for (int kb = bstart; kb < bend;
							             ++kb) {
								        int  jb  = b_cols_dev(kb);
								        auto bij = b_vals_dev(kb);
								        if constexpr (is_complex)
									        if (isConjTransB)
										        bij = Kokkos::
										            conj(
										                bij);
								        auto cij = aij * bij;
								        int  iy  = jb;
								        int  jy
								            = doTransA ? ia1 : ja;
								        acc += cij
								            * yin_dev(iy, jy);
							        }
						        }
					        } else {
						        // doTransB true: ix is jb (column index)
						        int jb     = ix;
						        int bstart = colPtrB_dev(jb);
						        int bend   = colPtrB_dev(jb + 1);
						        int len    = bend - bstart;
						        if (len > 32) {
							        Kokkos::parallel_for(
							            Kokkos::ThreadVectorRange(
							                member, len),
							            [&](const int t)
							            {
								            int pos = bstart + t;
								            int kpos
								                = idxListB_dev(pos);
								            int ib1 = b_row_of_dev(
								                kpos);
								            auto bij
								                = b_vals_dev(kpos);
								            if constexpr (
								                is_complex)
									            if (isConjTransB)
										            bij = Kokkos::
										                conj(
										                    bij);
								            auto cij = aij * bij;
								            int  iy  = ib1;
								            int  jy = doTransA ? ia1
								                               : ja;
								            acc += cij
								                * yin_dev(iy, jy);
							            });
						        } else {
							        for (int pos = bstart; pos < bend;
							             ++pos) {
								        int kpos
								            = idxListB_dev(pos);
								        int ib1
								            = b_row_of_dev(kpos);
								        auto bij = b_vals_dev(kpos);
								        if constexpr (is_complex)
									        if (isConjTransB)
										        bij = Kokkos::
										            conj(
										                bij);
								        auto cij = aij * bij;
								        int  iy  = ib1;
								        int  jy
								            = doTransA ? ia1 : ja;
								        acc += cij
								            * yin_dev(iy, jy);
							        }
						        }
					        }
				        }

				        results_dev((size_t)jx * (size_t)nrow_X + (size_t)ix) = acc;
			        });
		    });

		Exec().fence();

		// copy back results_dev to host and then into xout
		Kokkos::deep_copy(results_h, results_dev);
		for (int jx = 0; jx < nTarget; ++jx)
			for (int ix = 0; ix < nrow_X; ++ix)
				xout(ix, jx) = results_h[(size_t)jx * (size_t)nrow_X + (size_t)ix];

		Exec().fence();
#else
		int ia = 0;
		int ka = 0;
		int ib = 0;
		int kb = 0;
		for (ia = 0; ia < nrow_A; ia++) {
			int istarta = a.getRowPtr(ia);
			int ienda   = a.getRowPtr(ia + 1);
			for (ka = istarta; ka < ienda; ka++) {
				int               ja  = a.getCol(ka);
				ComplexOrRealType aij = a.getValue(ka);
				if (is_complex && isConjTransA) {
					aij = PsimagLite::conj(aij);
				};

				for (ib = 0; ib < nrow_B; ib++) {
					int istartb = b.getRowPtr(ib);
					int iendb   = b.getRowPtr(ib + 1);

					for (kb = istartb; kb < iendb; kb++) {
						int               jb  = b.getCol(kb);
						ComplexOrRealType bij = b.getValue(kb);
						if (is_complex && isConjTransB) {
							bij = PsimagLite::conj(bij);
						};

						ComplexOrRealType cij = aij * bij;

						int ix = (isTransB || isConjTransB) ? jb : ib;
						int jx = (isTransA || isConjTransA) ? ja : ia;
						int iy = (isTransB || isConjTransB) ? ib : jb;
						int jy = (isTransA || isConjTransA) ? ia : ja;

						xout(ix, jx) += cij * yin(iy, jy);
					};
				};
			};
		};
#endif
	}
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
