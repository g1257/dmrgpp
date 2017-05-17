#include "util.h"

template<typename ComplexOrRealType>
void den_kron_submatrix(const int nrow_A,
                        const int ncol_A,
                        const PsimagLite::Matrix<ComplexOrRealType>& a_,
                        const int nrow_B,
                        const int ncol_B,
                        const PsimagLite::Matrix<ComplexOrRealType>& b_,
                        const int nrindex,
                        const int ncindex,
                        const PsimagLite::Vector<int>::Type& rindex,
                        const PsimagLite::Vector<int>::Type& cindex,
                        PsimagLite::Matrix<ComplexOrRealType>& e_ )
{
/*
 * -------------------------------------------------
 * extract a submatrix out of kronecker product
 * equivalent to
 * C = kron(A,B), then E = C( rindex(:), cindex(:) )
 *
 * assume A, B are in dense matrix format
 * -------------------------------------------------
 */
	const int nrow_E = nrindex;
	const int ncol_E = ncindex;

	int ie = 0;
	int je = 0;

	/*
  * ------------------------
  * check rindex[], cindex[]
  * ------------------------
  */
#ifndef NDEBUG

	const int nrow_C = nrow_A * nrow_B;
	const int ncol_C = ncol_A * ncol_B;

	for(ie=0; ie < nrindex; ie++) {
		int ic = rindex[ie];
		int isok = (0 <= ic) && (ic < nrow_C);
		assert( isok );
	};

	for(je=0; je < ncindex; je++) {
		int jc = cindex[je];
		int isok = (0 <= jc) && (jc < ncol_C);
		assert( isok );
	};
#endif


	/*
  * --------------------------------
  * fill entries in E(ie,je) matrix
  * --------------------------------
  */

	for(je=0; je < ncol_E; je++) {
		for(ie=0; ie < nrow_E; ie++) {
			int ic = rindex[ie];
			int jc = cindex[je];
			/*
	 * --------------------------
	 * ic = [ib,ia] or
	 * ic =  ib + ia * nrow_B
	 *
	 * jc = [jb,ja] or
	 * jc =  jb + ja * ncol_B
	 * --------------------------
	 */
			int ib = MOD(ic,nrow_B);
			int ia = (ic - ib)/nrow_B;


			int jb = MOD(jc,ncol_B);
			int ja = (jc - jb)/ncol_B;



			assert((0 <= ia) && (ia < nrow_A) &&
			       (0 <= ib) && (ib < nrow_B) &&
			       ((ib + ia*nrow_B) == ic ));


			assert((0 <= ja) && (ja < ncol_A) &&
			       (0 <= jb) && (jb < ncol_B) &&
			       ((jb + ja*ncol_B) == jc));

			ComplexOrRealType cij = a_(ia,ja)*b_(ib,jb);

			e_(ie,je) = cij;
		};
	};
}
#undef A
#undef B
#undef E
