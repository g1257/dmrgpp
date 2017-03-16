#include "util.h"

void csr_matmul_pre( char trans_A, 
                     const int nrow_A,
                     const int ncol_A, 
                     const PsimagLite::Vector<int>::Type& arowptr,
                     const PsimagLite::Vector<int>::Type& acol,
                     const PsimagLite::Vector<double>::Type& aval,

                     const int nrow_Y, 
                     const int ncol_Y, 
                     const PsimagLite::MatrixNonOwned<const double>& yin,

                     const int nrow_X, 
                     const int ncol_X, 
                     PsimagLite::MatrixNonOwned<double>& xout)
{
   const int idebug = 0;
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
   PsimagLite::Matrix<double> a_(nrow_A, ncol_A);

   if (idebug >= 1) {
     printf("csr_matmul_pre:nrow_A %d,ncol_A %d,nnz_A %d\n",
             nrow_A, ncol_A, nnz_A );
     };

   csr2den( nrow_A, ncol_A,      arowptr, acol, aval,
            a_);
   
   den_matmul_pre( trans_A,
                   nrow_A, ncol_A,  a_,

                   nrow_Y,ncol_Y, yin,
                   nrow_X,ncol_X, xout);
   return;
   };


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

   int ia = 0;
   for(ia=0; ia < nrow_A; ia++) {
       int istart = arowptr[ia];
       int iend = arowptr[ia+1]-1;
       int k = 0;
       for(k=istart; k <= iend; k++) {
          int ja = acol[k];
          double aij = aval[k];
          double atji = aij;
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
          
   int ia = 0;
   for(ia=0; ia < nrow_A; ia++) {
       int istart = arowptr[ia];
       int iend = arowptr[ia+1]-1;
       int k = 0;
       for(k=istart; k <= iend; k++) {
          int ja = acol[k];
          double aij = aval[k];
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
#undef A
          




   

                   
