#include "util.h"

void csc_matmul_post( char trans_A, 
                     const int nrow_A,
                     const int ncol_A, 
                     const int acolptr[], 
                     const int arow[], 
                     const double aval[],

                     const int nrow_Y, 
                     const int ncol_Y, 
                     const double yin[],

                     const int nrow_X, 
                     const int ncol_X, 
                     double xout[] )

#define Y(iy,jy) yin[ (iy) + (jy)*nrow_Y ]
#define X(ix,jx) xout[ (ix) + (jx)*nrow_X ]

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

 int isTranspose = (trans_A == 'T') || (trans_A == 't');



 if (isTranspose) {
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
          double aij = aval[k];

          int ia = arow[k];
          assert((0 <= ia) && (ia < nrow_A));

          double atji = aij;
          
          int iy = 0;
          for(iy=0; iy < nrow_Y; iy++) {
            int ix = iy;
            int jx = ia;

            X(ix,jx) += (Y(iy,ja) * atji);
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
          double aij = aval[k];

          int ia = arow[k];
          assert((0 <= ia) && (ia < nrow_A));

          int iy = 0;
          for(iy = 0; iy < nrow_Y; iy++) {
              int ix = iy;
              int jx = ja;

              X(ix,jx) += ( Y(iy,ia) * aij );
              };
          };
       };
   }
}
          
#undef X
#undef Y
