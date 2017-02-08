#include "util.h"

void csc_matmul_pre( char trans_A, 
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

#define Y(iy,jy)  yin[ (iy) + (jy)*nrow_Y ]
#define X(ix,jx) xout[ (ix) + (jx)*nrow_X ]

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
         
          double aij = aval[k];
          double atji = aij;
          int jy = 0;
          for(jy=0; jy < ncol_Y; jy++) {
            int ix = ja;
            int jx = jy;
            X(ix,jx) += (atji * Y(ia,jy));
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
          double aij = aval[k];

          int ia = arow[k];
          assert((0 <= ia) && (ia < nrow_A));

          int jy = 0;
          for(jy=0; jy < ncol_Y; jy++) {
            int ix = ia;
            int jx = jy;

            X(ix,jx) += (aij * Y(ja,jy));
            };
         };
      };
  };
}
          
#undef X
#undef Y
