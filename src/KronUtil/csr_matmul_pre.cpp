#include "util.h"

template<typename ComplexOrRealType>
void csr_matmul_pre( char trans_A, 
                     const PsimagLite::CrsMatrix<ComplexOrRealType>& a,
                     const int nrow_Y,
                     const int ncol_Y,
                     const PsimagLite::MatrixNonOwned<const ComplexOrRealType>& yin,
                     const int nrow_X,
                     const int ncol_X,
                     PsimagLite::MatrixNonOwned<ComplexOrRealType>& xout)
{
/*
 * -------------------------------------------------------
 * A in compressed sparse ROW format
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

	const int nrow_A = a.rows();
	int isTranspose = (trans_A == 'T') || (trans_A == 't');

	if (isTranspose) {
	/*
	 *   ----------------------------------------------------------
	 *   X(nrow_X,ncol_X) +=  tranpose(A(nrow_A,ncol_A))*Y(nrow_Y,ncol_Y)
	 *   X(ix,jx) +=  transpose( A(ia,ja) ) * Y(iy,jy)
	 *   X(ja,jy) +=  sum( At(ja, ia) * Y(ia,jy), over ia)
	 *   ----------------------------------------------------------
	 */

		assert(static_cast<SizeType>(nrow_X) == a.cols());
		assert(nrow_A == nrow_Y && ncol_X == ncol_Y);

		int ia = 0;
		for(ia=0; ia < nrow_A; ia++) {
			int istart = a.getRowPtr(ia);
			int iend = a.getRowPtr(ia + 1);
			int k = 0;
			for(k=istart; k < iend; k++) {
				int ja = a.getCol(k);
				ComplexOrRealType aij = a.getValue(k);
				ComplexOrRealType atji = aij;
				int jy = 0;
				for(jy=0; jy < ncol_Y; jy++) {
					int ix = ja;
					int jx = jy;
					xout(ix,jx) += (atji * yin(ia,jy));
				}
			}
		}
	} else  {
	/*
	 * ---------------------------------------------
	 * X(nrow_X,ncol_X) += A(nrow_A,ncol_A) * Y(nrow_Y,ncol_Y)
	 * X(ia,jy) += sum( A(ia,ja)*Y(ja,jy), over ja )
	 * ---------------------------------------------
	 */
		assert(nrow_X == nrow_A);
		assert(a.cols() == static_cast<SizeType>(nrow_Y) && (ncol_X == ncol_Y));

		int ia = 0;
		for(ia=0; ia < nrow_A; ia++) {
			int istart = a.getRowPtr(ia);
			int iend = a.getRowPtr(ia + 1);
			int k = 0;
			for(k=istart; k < iend; k++) {
				int ja = a.getCol(k);
				ComplexOrRealType aij = a.getValue(k);
				int jy = 0;

				for(jy=0; jy < ncol_Y; jy++) {
					int ix = ia;
					int jx = jy;

					xout(ix,jx) += (aij * yin(ja,jy));
				}
			}
		}
	}
}









