#include "util.h"

template<typename ComplexOrRealType>
void csc_matmul_pre(char trans_A,
                    const int nrow_A,
                    const int ncol_A,
                    const PsimagLite::Vector<int>::Type& acolptr,
                    const PsimagLite::Vector<int>::Type& arow,
                    const typename PsimagLite::Vector<ComplexOrRealType>::Type& aval,
                    const int nrow_Y,
                    const int ncol_Y,
                    const PsimagLite::Matrix<ComplexOrRealType>& yin,
                    const int nrow_X,
                    const int ncol_X,
                    PsimagLite::Matrix<ComplexOrRealType>& xout )
{
	/*
 * -------------------------------------------------------
 * A in compressed sparse COLUMN format
 *
 * compute   X +=  op(A) * Y
 * where op(A) is transpose(A)   if trans_A = 'T' or 't'
 *       op(A) is A              otherwise
 *
 * if need transpose(A) then
 *   X(nrow_X,ncol_X) +=  tranpose(A(nrow_A,ncol_A))*Y(nrow_Y,ncol_Y)
 *   requires nrow_X == ncol_A, ncol_X == ncol_Y, nrow_A == nrow_Y
 *
 * if need A then
 *  X(nrow_X,ncol_X) += A(nrow_A,ncol_A) * Y(nrow_Y,ncol_Y)
 *  requires  nrow_X == nrow_A, ncol_A == nrow_Y, ncol_X == ncol_Y
 * -------------------------------------------------------
 */

	int isTranspose = (trans_A == 'T') || (trans_A == 't');



	if (isTranspose) {
		/*
	*   ----------------------------------------------------------
	*   X(nrow_X,ncol_X) +=  tranpose(A(nrow_A,ncol_A))*Y(nrow_Y,ncol_Y)
	*   X(ix,jx) +=  transpose( A(ia,ja) ) * Y(iy,jy)
	*   X(ja,jy) +=  sum( At(ja, ia) * Y(ia,jy), over ia)
	*   ----------------------------------------------------------
	*/

		assert((nrow_X == ncol_A) &&  (nrow_A == nrow_Y) && (ncol_X == ncol_Y));

		int ja = 0;
		for(ja=0; ja < ncol_A; ja++) {
			int istart = acolptr[ja];
			int iend = acolptr[ja+1]-1;
			int k = 0;
			for(k=istart; k <= iend; k++) {
				int ia = arow[k];
				assert((0 <= ia) && (ia < nrow_A));

				ComplexOrRealType aij = aval[k];
				ComplexOrRealType atji = aij;
				int jy = 0;
				for(jy=0; jy < ncol_Y; jy++) {
					int ix = ja;
					int jx = jy;
					xout(ix,jx) += (atji * yin(ia,jy));
				};
			};
		};
	}
	else  {
		/*
	* ---------------------------------------------
	* X(nrow_X,ncol_X) += A(nrow_A,ncol_A) * Y(nrow_Y,ncol_Y)
	* X(ia,jy) += sum( A(ia,ja)*Y(ja,jy), over ja )
	* ---------------------------------------------
	*/
		assert((nrow_X == nrow_A) && (ncol_A == nrow_Y) && (ncol_X == ncol_Y));

		int ja = 0;
		for(ja=0; ja < ncol_A; ja++) {
			int istart = acolptr[ja];
			int iend = acolptr[ja+1]-1;
			int k = 0;
			for(k=istart; k <= iend; k++) {
				ComplexOrRealType aij = aval[k];

				int ia = arow[k];
				assert((0 <= ia) && (ia < nrow_A));

				int jy = 0;
				for(jy=0; jy < ncol_Y; jy++) {
					int ix = ia;
					int jx = jy;

					xout(ix,jx) += (aij * yin(ja,jy));
				};
			};
		};
	};
}

#undef X
#undef Y
