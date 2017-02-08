#include "util.h"

void csr_matmul_post( char trans_A, 
                     const int nrow_A,
                     const int ncol_A, 
                     const int arowptr[], 
                     const int acol[], 
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
    const int idebug = 0;
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


 const double threshold = 0.5;

 int nnz_A = csr_nnz( nrow_A, arowptr );
 int use_dense_storage = (nnz_A >= threshold * nrow_A * ncol_A);

 if (use_dense_storage) {
   /*
    * ----------------------------------
    * sparse matrix A is sufficiently dense that
    * we convert to dense storage to perform 
    * more efficient matrix operations
    * ----------------------------------
    */
   double a_[nrow_A * ncol_A];
#define A(ia,ja)  a_[ (ia) + (ja)*nrow_A ]
   if (idebug >= 1) {
     printf("csr_matmul_post:nrow_A %d,ncol_A %d,nnz_A %d\n",
             nrow_A, ncol_A, nnz_A );
     };

   csr2den( nrow_A, ncol_A,         arowptr, acol, aval,
            &(A(0,0))  );
   
   den_matmul_post( trans_A,
                   nrow_A, ncol_A,  &(A(0,0)),

                   nrow_Y,ncol_Y, &(Y(0,0)),
                   nrow_X,ncol_X, &(X(0,0)) );

   return;
   };

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

   int ia = 0;
   for(ia=0; ia < nrow_A; ia++) {
       int istart = arowptr[ia];
       int iend = arowptr[ia+1]-1;
       int k = 0;
       for(k=istart; k <= iend; k++) {
          int ja = acol[k];
          double aij = aval[k];
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
          
   int ia = 0;
   for(ia=0; ia < nrow_A; ia++) {
       int istart = arowptr[ia];
       int iend = arowptr[ia+1]-1;
       int k = 0;
       for(k=istart; k <= iend; k++) {
          int ja = acol[k];
          double aij = aval[k];
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
