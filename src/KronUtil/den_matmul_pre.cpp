#include "util.h"

void den_matmul_pre( const char trans_A, 
                     const int nrow_A,
                     const int ncol_A, 
                     const PsimagLite::Matrix<double>& a_,

                     const int nrow_Y, 
                     const int ncol_Y, 
                     const PsimagLite::Matrix<double>& yin,

                     const int nrow_X, 
                     const int ncol_X, 
                     PsimagLite::Matrix<double>& xout)

{
/*
 * -------------------------------------------------------
 * A in dense matrix format
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

 const int use_blas =  TRUE;
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

    if (use_blas) {
        /*
         * ---------------------
         * X  = transpose(A)*Y + X
         * ---------------------
         */
         char trans1 = 'T';
         char trans2 = 'N';
         int mm = nrow_X;
         int nn = ncol_X;
         int kk = nrow_Y;

         double alpha = 1;
         double beta = 1;
         int ld1 = nrow_A;
         int ld2 = nrow_Y;
         int ld3 = nrow_X;

         dgemm_( &trans1, &trans2,
                 &mm, &nn, &kk,
                 &alpha, &(a_(0,0)), &ld1,
                         &(yin(0,0)), &ld2,
                 &beta,  &(xout(0,0)), &ld3 );


      }
    else {
      int jx = 0;
      
      for(jx=0; jx < ncol_X; jx++) {
       int ix = 0;
       for(ix=0; ix < nrow_X; ix++) {
          int ja = ix;
          int jy = jx;

          double dsum = 0;
          int ia = 0;
          for(ia=0; ia < nrow_A; ia++) {
              int iy = ia;
              double aij = a_(ia,ja);
              double atji = aij;
              dsum +=  (atji * yin(iy,jy));
              };
          xout(ix,jx) += dsum;
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

    if (use_blas) {
       /*
        * --------------
        * X =  A * Y + X
        * --------------
        */
        char trans1 = 'N';
        char trans2 = 'N';
        int mm = nrow_X;
        int nn = ncol_X;
        int kk = nrow_Y;

        double alpha = 1;
        double beta = 1;
        int ld1 = nrow_A;
        int ld2 = nrow_Y;
        int ld3 = nrow_X;

        dgemm_( &trans1, &trans2,
                &mm, &nn, &kk,
                &alpha,  &(a_(0,0)), &ld1,
                         &(yin(0,0)), &ld2,
                &beta,   &(xout(0,0)), &ld3 );
       }
    else {
       /*
        * ----------------------------------------------
        * X(ix,jx) += sum( A(ia,ja) * Y(ja,jy), over ja)
        * ----------------------------------------------
        */
       int jx = 0;
       for(jx=0; jx < ncol_X; jx++) {
         int ix = 0;
         for(ix=0; ix < nrow_X; ix++) {
           int ia = ix;
           int jy = jx;

           double dsum = 0;
           int ja = 0;
           for(ja=0; ja < ncol_A; ja++) {
              int iy = ja;
              double aij = a_(ia,ja);
              dsum += (aij * yin(iy,jy));
              };
           xout(ix,jx) += dsum;
           };
         };
              
    };
          
  };
}
          
#undef X
#undef Y
#undef A
