#include "util.h"

template<typename ComplexOrRealType>
void den_matmul_post(const char trans_A,
                     const int nrow_A,
                     const int ncol_A,
                     const PsimagLite::Matrix<ComplexOrRealType>& a_,
                     const int nrow_Y,
                     const int ncol_Y,
                     const PsimagLite::MatrixNonOwned<const ComplexOrRealType>& yin,
                     const int nrow_X,
                     const int ncol_X,
                     PsimagLite::MatrixNonOwned<ComplexOrRealType>& xout)
{
	const bool is_complex = PsimagLite::IsComplexNumber<ComplexOrRealType>::True;
/*
 * -------------------------------------------------------
 * A in dense matrix format
 *
 * compute   X +=  Y * op(A)
 * where op(A) is transpose(A)   if trans_A = 'T' or 't'
 *       op(A) is conj(transpose(A)) if trans_A = 'C' or 'c'
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

	int isTranspose = (trans_A == 'T') || (trans_A == 't');
	int isConjTranspose = (trans_A == 'C') || (trans_A == 'c');

	const bool use_blas = true;

	if (isTranspose || isConjTranspose) {
	/*
	*   ----------------------------------------------------------
	*   X(nrow_X,ncol_X) +=  Y(nrow_Y,ncol_Y) * transpose(A(nrow_A,ncol_A))
	*   X(ix,jx) +=  Y(iy,jy) * transpose( A(ia,ja) )
	*   X(ix,jx) += sum( Y(iy, ja) * At(ja,ia), over ja )
	*   ----------------------------------------------------------
	*/

		assert((nrow_X == nrow_Y) && (ncol_Y == ncol_A) && (ncol_X == nrow_A));

		if (use_blas) {
		/*
		 * ------------------------
		 * X = Y * transpose(A) + X
		 * X(ix,jx) += sum( Y(ix,ja) * A(ia,ja), over ja)
		 * ------------------------
		 */
			char trans1 = 'N';
			char trans2 = trans_A;



			int mm = nrow_X;
			int nn = ncol_X;
			int kk = ncol_Y;
			ComplexOrRealType alpha = 1;
			ComplexOrRealType beta = 1;
			int ld1 = nrow_Y;
			int ld2 = nrow_A;
			int ld3 = nrow_X;

			psimag::BLAS::GEMM(trans1, trans2,
			                      mm, nn, kk,
			                      alpha,  &(yin(0,0)), ld1,
			                      &(a_(0,0)), ld2,
			                      beta,   &(xout(0,0)), ld3);

		}
		else {

			int jx = 0;
			for(jx=0; jx < ncol_X; jx++) {
				int ix = 0;
				for(ix=0; ix < nrow_X; ix++) {
					ComplexOrRealType dsum = 0;
					int iy = ix;
					int ja = 0;
					for(ja=0; ja < ncol_A; ja++) {
						int jy = ja;
						int ia = jx;
						ComplexOrRealType aij = a_(ia,ja);
						ComplexOrRealType atji = aij;
						if (is_complex && isConjTranspose) {
							atji = PsimagLite::conj( atji );
						};

						dsum += (yin(iy,jy) * atji);
					};
					xout(ix,jx) += dsum;
				};
			};


		};

	}
	else  {
	/*
	* ---------------------------------------------
	* X(nrow_X,ncol_X) += Y(nrow_Y,ncol_Y) * A(nrow_A,ncol_A)
	* X(ix,jx) += sum( Y(iy,ia) * A(ia,ja), over ia )
	* ---------------------------------------------
	*/
		assert((nrow_X == nrow_Y) && (ncol_Y == nrow_A) && (ncol_X == ncol_A));

		if (use_blas) {
		/*
		 * -------------
		 * X = Y * A + X
		 * -------------
		 */
			char trans1 = 'N';
			char trans2 = 'N';
			int mm = nrow_X;
			int nn = ncol_X;
			int kk = ncol_Y;
			ComplexOrRealType alpha = 1;
			ComplexOrRealType beta = 1;
			int ld1 = nrow_Y;
			int ld2 = nrow_A;
			int ld3 = nrow_X;

			psimag::BLAS::GEMM(trans1, trans2,
			                      mm, nn, kk,
			                      alpha,  &(yin(0,0)), ld1,
			                      &(a_(0,0)), ld2,
			                      beta,   &(xout(0,0)), ld3);

		}
		else {
			int jx = 0;
			for(jx=0; jx < ncol_X; jx++) {
				int ix = 0;
				for(ix=0; ix < ncol_X; ix++) {
					ComplexOrRealType dsum = 0;
					int ia = 0;
					for(ia=0; ia < nrow_A; ia++) {
						int ja = jx;
						int iy = ix;
						int jy = ia;
						ComplexOrRealType aij = a_(ia,ja);

						dsum += (aij * yin(iy,jy));
					};
					xout(ix,jx) += dsum;
				};
			};
		};


	}
}

#undef X
#undef Y
#undef A
