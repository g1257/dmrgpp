#include "util.h"

template<typename ComplexOrRealType>
void csc_matmul_post(char trans_A,
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
                     PsimagLite::Matrix<ComplexOrRealType>& xout)
{
/*
 * -------------------------------------------------------
 * A in compress sparse COLUMN format
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
 const bool is_complex = std::is_same<ComplexOrRealType,std::complex<double> >::value ||
	                 std::is_same<ComplexOrRealType,std::complex<float> >::value;

 int isTranspose = (trans_A == 'T') || (trans_A == 't');
 int isConjTranspose = (trans_A == 'C') || (trans_A == 'c');



 if (isTranspose || isConjTranspose) {
   /*
    *   ----------------------------------------------------------
    *   X(nrow_X,ncol_X) +=  Y(nrow_Y,ncol_Y) * transpose(A(nrow_A,ncol_A))
    *   X(ix,jx) +=  Y(iy,jy) * transpose( A(ia,ja) )
    *   X(ix,jx) += sum( Y(iy, ja) * At(ja,ia), over ja )
    *   ----------------------------------------------------------
    */

    assert((nrow_X == nrow_Y) && (ncol_Y == ncol_A) && (ncol_X == nrow_A));

   int ja = 0;
   for(ja=0; ja < ncol_A; ja++) {
       int istart = acolptr[ja];
       int iend = acolptr[ja+1]-1;
       int k = 0;
       for(k=istart; k <= iend; k++) {
          ComplexOrRealType aij = aval[k];

          int ia = arow[k];
          assert((0 <= ia) && (ia < nrow_A));

          ComplexOrRealType atji =  aij;
	  if (is_complex && isConjTranspose) {
		  atji = PsimagLite::conj(atji);
	  };

          int iy = 0;
          for(iy=0; iy < nrow_Y; iy++) {
            int ix = iy;
            int jx = ia;

            xout(ix,jx) += (yin(iy,ja) * atji);
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

   int ja = 0;
   for(ja=0; ja < ncol_A; ja++) {
       int istart = acolptr[ja];
       int iend = acolptr[ja+1]-1;
       int k = 0;
       for(k=istart; k <= iend; k++) {
          ComplexOrRealType aij = aval[k];

          int ia = arow[k];
          assert((0 <= ia) && (ia < nrow_A));

          int iy = 0;
          for(iy = 0; iy < nrow_Y; iy++) {
              int ix = iy;
              int jx = ja;

              xout(ix,jx) += (yin(iy,ia) * aij );
              };
          };
       };
   }
}

#undef X
#undef Y
