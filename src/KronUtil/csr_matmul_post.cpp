#include "util.h"

template <typename ComplexOrRealType>
void csr_matmul_post(char trans_A,
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
	const bool is_complex = PsimagLite::IsComplexNumber<ComplexOrRealType>::True;
	const int nrow_A = a.rows();
	int isTranspose = (trans_A == 'T') || (trans_A == 't');
	int isConjTranspose = (trans_A == 'C') || (trans_A == 'c');
	int isConj = (trans_A == 'Z') || (trans_A == 'z');

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

		int ia = 0;
		for (ia = 0; ia < nrow_A; ia++) {
			int istart = a.getRowPtr(ia);
			int iend = a.getRowPtr(ia + 1);
			int k = 0;
			for (k = istart; k < iend; k++) {
				int ja = a.getCol(k);
				ComplexOrRealType aij = a.getValue(k);
				ComplexOrRealType atji = aij;
				if (is_complex && isConjTranspose) {
					atji = PsimagLite::conj(atji);
				};

				int iy = 0;
				for (iy = 0; iy < nrow_Y; iy++) {
					int ix = iy;
					int jx = ia;

					xout(ix, jx) += (yin(iy, ja) * atji);
				}
			}
		}
	} else {
		/*
		 * ---------------------------------------------
		 * X(nrow_X,ncol_X) += Y(nrow_Y,ncol_Y) * A(nrow_A,ncol_A)
		 * X(ix,jx) += sum( Y(iy,ia) * A(ia,ja), over ia )
		 * ---------------------------------------------
		 */
		assert(nrow_X == nrow_Y);
		assert(ncol_Y == nrow_A && static_cast<SizeType>(ncol_X) == a.cols());

		int ia = 0;
		for (ia = 0; ia < nrow_A; ia++) {
			int istart = a.getRowPtr(ia);
			int iend = a.getRowPtr(ia + 1);
			int k = 0;
			for (k = istart; k < iend; k++) {
				int ja = a.getCol(k);
				ComplexOrRealType aij = a.getValue(k);
				if (is_complex && isConj) {
					aij = PsimagLite::conj(aij);
				};

				int iy = 0;
				for (iy = 0; iy < nrow_Y; iy++) {
					int ix = iy;
					int jx = ja;

					xout(ix, jx) += (yin(iy, ia) * aij);
				}
			}
		}
	}
}
