#include "util.h"
#include <PsimagLite/KokkosType.h>

#include <Kokkos_Core.hpp>
#include <Kokkos_Profiling_ScopedRegion.hpp>

#include <Kokkos_Profiling_ScopedRegion.hpp>

#include <Kokkos_Profiling_ScopedRegion.hpp>

template <typename ComplexOrRealType>
void csr_matmul_post(char                                                       trans_A,
                     const PsimagLite::CrsMatrix<ComplexOrRealType>&            a,
                     const int                                                  nrow_Y,
                     const int                                                  ncol_Y,
                     const PsimagLite::MatrixNonOwned<const ComplexOrRealType>& yin,
                     const int                                                  nrow_X,
                     const int                                                  ncol_X,
                     PsimagLite::MatrixNonOwned<ComplexOrRealType>&             xout)
{
	Kokkos::Profiling::ScopedRegion region("PsimagLite::csr_matmul_post");
	/*
	 * -------------------------------------------------------
	 * A in compressed sparse ROW format
	 *
	 * compute   X +=  Y * op(A)
	 * where op(A) is transpose(A)   if trans_A = 'T' or 't'
	 *       op(A) is A              otherwise
	 *
	 * if need transpose(A) then
	 *   X(nrow_X,ncol_X) +=  Y(nrow_Y,ncol_Y) * tranpose(A(nrow_A,ncol_A))
	 *   requires (nrow_X == nrow_Y) && (ncol_Y == ncol_A) && (ncol_X == nrow_A)
	 *
	 * if need A then
	 *  X(nrow_X,ncol_X) +=  Y(nrow_Y,ncol_Y) * A(nrow_A,ncol_A)
	 *  requires  (nrow_X == nrow_Y) && ( ncol_Y == nrow_A) && (ncol_X == ncol_A)
	 * -------------------------------------------------------
	 */
	const bool is_complex      = PsimagLite::IsComplexNumber<ComplexOrRealType>::True;
	const int  nrow_A          = a.rows();
	int        isTranspose     = (trans_A == 'T') || (trans_A == 't');
	int        isConjTranspose = (trans_A == 'C') || (trans_A == 'c');
	int        isConj          = (trans_A == 'Z') || (trans_A == 'z');

	using ExecutionSpace = Kokkos::DefaultExecutionSpace;
	ExecutionSpace exec;

	using KokkosScalar = typename PsimagLite::KokkosType<ComplexOrRealType>::type;

	const int nnz = a.nonZeros();

	Kokkos::View<const int*, Kokkos::HostSpace, Kokkos::MemoryUnmanaged> rowptr_host(
	    &a.getRowPtr(0), nrow_A + 1);
	Kokkos::View<const int*, Kokkos::HostSpace, Kokkos::MemoryUnmanaged> cols_host(&a.getCol(0),
	                                                                               nnz);
	Kokkos::View<const KokkosScalar*, Kokkos::HostSpace, Kokkos::MemoryUnmanaged> vals_host(
	    reinterpret_cast<const KokkosScalar*>(&a.getValue(0)), nnz);

	Kokkos::View<KokkosScalar**> x_dev_out("x_dev_out", nrow_Y, ncol_X);
	Kokkos::View<const KokkosScalar**,
	             Kokkos::LayoutLeft,
	             Kokkos::HostSpace,
	             Kokkos::MemoryUnmanaged>
	     yin_host(reinterpret_cast<const KokkosScalar*>(&yin(0, 0)), nrow_Y, ncol_Y);
	auto y_dev = Kokkos::create_mirror_view_and_copy(ExecutionSpace {}, yin_host);

	// Copy CSR arrays to device so the TeamPolicy kernel can access them
	auto d_rowptr = Kokkos::create_mirror_view_and_copy(ExecutionSpace {}, rowptr_host);
	auto d_cols   = Kokkos::create_mirror_view_and_copy(ExecutionSpace {}, cols_host);
	auto d_vals   = Kokkos::create_mirror_view_and_copy(ExecutionSpace {}, vals_host);

	// FIXME We can try using KokkosKernels directly later
	using team_policy = Kokkos::TeamPolicy<ExecutionSpace>;
	using member_type = team_policy::member_type;

	// One team per output row iy; let Kokkos pick team/vector sizes
	team_policy policy(nrow_Y, Kokkos::AUTO);
	if (isTranspose || isConjTranspose) {
		/*
		 *   ----------------------------------------------------------
		 *   X(nrow_X,ncol_X) +=  Y(nrow_Y,ncol_Y) * transpose(A(nrow_A,ncol_A))
		 *   X(ix,jx) +=  Y(iy,jy) * transpose( A(ia,ja) )
		 *   X(ix,jx) += sum( Y(iy, ja) * At(ja,ia), over ja )
		 *   ----------------------------------------------------------
		 */

		assert(nrow_X == nrow_Y);
		assert(static_cast<SizeType>(ncol_Y) == a.cols() && (ncol_X == nrow_A));

		Kokkos::parallel_for(
		    "csr_matmul_post::team_transpose",
		    policy,
		    KOKKOS_LAMBDA(const member_type& team) {
			    const int iy = team.league_rank();

			    // create a subview for the current Y row for faster access
			    auto yrow = Kokkos::subview(y_dev, iy, Kokkos::ALL);

			    // For transpose: X(iy,ia) += sum_j Y(iy,ja)*A(ia,ja)
			    Kokkos::parallel_for(
			        Kokkos::TeamThreadRange(team, nrow_A),
			        [&](int ia)
			        {
				        int          istart = d_rowptr(ia);
				        int          iend   = d_rowptr(ia + 1);
				        KokkosScalar local_sum;
				        // reduce contributions across the nonzeros of row
				        // ia
				        Kokkos::parallel_reduce(
				            Kokkos::ThreadVectorRange(team, istart, iend),
				            [&](int k, KokkosScalar& lsum)
				            {
					            int          ja  = d_cols(k);
					            KokkosScalar aij = d_vals(k);
					            if constexpr (is_complex)
						            if (isConjTranspose)
							            aij = Kokkos::conj(aij);
					            lsum += yrow(ja) * aij;
				            },
				            local_sum);
				        Kokkos::single(Kokkos::PerThread(team),
				                       [&]() { x_dev_out(iy, ia) += local_sum; });
			        });
		    });
	} else {
		/*
		 * ---------------------------------------------
		 * X(nrow_X,ncol_X) += Y(nrow_Y,ncol_Y) * A(nrow_A,ncol_A)
		 * X(ix,jx) += sum( Y(iy,ia) * A(ia,ja), over ia )
		 * ---------------------------------------------
		 */
		assert(nrow_X == nrow_Y);
		assert(ncol_Y == nrow_A && static_cast<SizeType>(ncol_X) == a.cols());

		Kokkos::parallel_for(
		    "csr_matmul_post::team_no_transpose",
		    policy,
		    KOKKOS_LAMBDA(const member_type& team) {
			    const int iy = team.league_rank();

			    // create a subview for the current Y row for faster access
			    auto yrow = Kokkos::subview(y_dev, iy, Kokkos::ALL);

			    // Non-transpose: X(iy,ja) += sum_ia Y(iy,ia)*A(ia,ja)
			    // Parallelize over ia (TeamThreadRange), then vectorize over
			    // nonzeros and use atomics for updates
			    Kokkos::parallel_for(
			        Kokkos::TeamThreadRange(team, nrow_A),
			        [&](int ia)
			        {
				        int          istart = d_rowptr(ia);
				        int          iend   = d_rowptr(ia + 1);
				        KokkosScalar yval   = yrow(ia);
				        Kokkos::parallel_for(
				            Kokkos::ThreadVectorRange(team, istart, iend),
				            [&](int k)
				            {
					            int          ja  = d_cols(k);
					            KokkosScalar aij = d_vals(k);
					            if constexpr (is_complex)
						            if (isConj)
							            aij = Kokkos::conj(aij);
					            KokkosScalar prod = yval * aij;
					            Kokkos::atomic_add(&x_dev_out(iy, ja), prod);
				            });
			        });
		    });
	}
	auto xhost = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace {}, x_dev_out);
	for (int iy = 0; iy < nrow_Y; ++iy) {
		for (int jx = 0; jx < ncol_X; ++jx)
			xout(iy, jx) += static_cast<ComplexOrRealType>(xhost(iy, jx));
	}
}
