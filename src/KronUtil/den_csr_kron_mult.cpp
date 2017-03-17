#include "util.h"

void den_csr_kron_mult_method(const int imethod,
                    const char transA,
                    const char transB,

                    const PsimagLite::Matrix<double>& a_,

                    const PsimagLite::CrsMatrix<double>& b,

                    const double* yin_,
                          double* xout_)
{
     const int isTransA = (transA == 'T') || (transA == 't');
     const int isTransB = (transB == 'T') || (transB == 't');
    
	 const int nrow_A = a_.n_row();
	 const int ncol_A = a_.n_col();
	 const int nrow_B = b.row();
	 const int ncol_B = b.col();

     const int nrow_1 = (isTransA) ? ncol_A : nrow_A;
     const int ncol_1 = (isTransA) ? nrow_A : ncol_A;
     const int nrow_2 = (isTransB) ? ncol_B : nrow_B;
     const int ncol_2 = (isTransB) ? nrow_B : ncol_B;

     const int nrow_X = nrow_2;
     const int ncol_X = nrow_1;
     const int nrow_Y = ncol_2;
     const int ncol_Y = ncol_1;

	 PsimagLite::MatrixNonOwned<const double> yin(nrow_Y, ncol_Y, yin_);
	 PsimagLite::MatrixNonOwned<double> xout(nrow_X, ncol_X, xout_);

     assert((imethod == 1) ||
                (imethod == 2) ||
                (imethod == 3));

     int nnz_A = den_nnz( nrow_A, ncol_A, a_);
     int nnz_B = csr_nnz(b);
     int has_work = (nnz_A >= 1) && (nnz_B >= 1);
     if (!has_work) {
         return;
         };
/*
 *   -------------------------------------------------------------
 *   A in dense matrix storage format
 *   B in compressed sparse ROW format
 *
 *   X += kron( op(A), op(B)) * Y 
 *   X +=  op(B) * Y * transpose(op(A))
 * 
 *   nrow_X = nrow_2,   ncol_X = nrow_1
 *   nrow_Y = ncol_2,   nrow_Y = nrow_2
 *
 *   that can be computed as either
 *   imethod == 1
 *
 *   X(ib,ia) +=  (B(ib,jb) * Y( jb,ja)) * transpose( A(ia,ja))
 *
 *   X(ix,jx) +=  (B2(ib2,jb2) * Y(iy,jy)) * transpose(A1(ia1,ja1))
 *
 *
 *
 *   X(ib,ia) += (B(ib,jb) * Y(jb,ja) ) * transpose(A(ia,ja)   or
 *               BY(ib,ja) = B(ib,jb)*Y(jb,ja)
 *               BY is nrow_B by ncol_A, need   2*nnz(B)*ncolA flops               
 *
 *   X(ib,ia) +=   BY(ib,ja) * transpose(A(ia,ja)) need 2*nnz(A)*nrowB flops
 *
 *   imethod == 2
 *
 *   X(ib,ia) += B(ib,jb) * (Y(jb,ja) * transpose(A))    or
 *                YAt(jb,ia) = Y(jb,ja) * transpose(A(ia,ja))    
 *                YAt is ncolB by nrowA, need 2*nnz(A) * ncolB flops
 *
 *   X(ib,ia) += B(ib,jb) * YAt(jb,ia)  need nnz(B) * nrowA flops
 *
 *   imethod == 3
 *
 *   X += kron(A,B) * Y   by visiting all non-zero entries in A, B        
 *
 *   this is feasible only if A and B are very sparse, need nnz(A)*nnz(B) flops
 *   -------------------------------------------------------------
 */



 if (imethod == 1) {

    /*
     *  --------------------------------------------
     *  BY(ib,ja) = (B(ib,jb))*Y(jb,ja)
     *
     *  X(ib,ia) += BY(ib,ja ) * transpose(A(ia,ja))
     *  --------------------------------------------
     */

    int nrow_BY = nrow_X;
    int ncol_BY = ncol_Y;
    PsimagLite::Matrix<double> by_(nrow_BY, ncol_BY);
	PsimagLite::MatrixNonOwned<double> byRef(by_);
	PsimagLite::MatrixNonOwned<const double> byConstRef(by_);

    /*
     * ---------------
     * setup BY(nrow_BY,ncol_BY)
     * ---------------
     */

    {
     int iby = 0;
     int jby = 0;

     // not needed FIXME
     for(jby=0; jby < ncol_BY; jby++) {
     for(iby=0; iby < nrow_BY; iby++) {
        by_(iby,jby) = 0;
        };
        };
     }



   {
    /*
     * ------------------------------
     * BY(ib,ja)  = B(ib,jb)*Y(jb,ja)
     * ------------------------------
     */
    const char trans = (isTransB) ? 'T' : 'N';
    csr_matmul_pre(  trans,
                     b,

                     nrow_Y,
                     ncol_Y,
                     yin,

                     nrow_BY,
                     ncol_BY,
                     byRef);

   }

   {
    /*
     * -------------------------------------------
     * X(ib,ia) += BY(ib,ja) * transpose(A(ia,ja))
     * -------------------------------------------
     */
     const char trans = (isTransA) ? 'N' : 'T';

     den_matmul_post( 
                     trans,
                     nrow_A,
                     ncol_A,
                     a_,

                     nrow_BY,
                     ncol_BY,
                     byConstRef,

                     nrow_X,
                     ncol_X,
                     xout);

    }
   }
 else if (imethod == 2) {
    /*
     * ---------------------
     * YAt(jb,ia) = Y(jb,ja) * tranpose(A(ia,ja))
     * X(ib,ia) += B(ib,jb) * YAt(jb,ia)
     * ---------------------
     */

    int nrow_YAt = nrow_Y;
    int ncol_YAt = ncol_X;
    PsimagLite::Matrix<double> yat_(nrow_YAt, ncol_YAt);
	PsimagLite::MatrixNonOwned<const double> yatConstRef(yat_);
	PsimagLite::MatrixNonOwned<double> yatRef(yat_);

    /*
     * ----------------
     * setup YAt(jb,ia)
     * ----------------
     */

    {
    int iy = 0;
    int jy = 0;
  
    for(jy=0; jy < ncol_YAt; jy++) {
    for(iy=0; iy < nrow_YAt; iy++) {
       yat_(iy,jy) = 0;
       };
       };
    }
    


   
    {
     /*
      * ---------------------
      * YAt(jb,ia) = Y(jb,ja) * tranpose(A(ia,ja)
      * ---------------------
      */
    const char transa = (isTransA) ? 'N' : 'T';



    den_matmul_post( transa,
                     nrow_A,
                     ncol_A,
                     a_,

                     nrow_Y, 
                     ncol_Y,
                     yin,
    
                     nrow_YAt, 
                     ncol_YAt,
                     yatRef);
     }




    {
    /*
     * ------------
     * X(ib,ia) += B(ib,jb) * YAt(jb,ia)
     * ------------
     */

    const char trans = (isTransB) ? 'T' : 'N';
    csr_matmul_pre( trans,
                    b,

                    nrow_YAt,
                    ncol_YAt,
                    yatConstRef,

                     nrow_X,
                     ncol_X, 
                     xout);
                     
      }
   }
 else if (imethod == 3) {
   /*
    * ---------------------------------------------
    * C = kron(A,B)
    * C([ib,ia], [jb,ja]) = A(ia,ja)*B(ib,jb)
    * X([ib,ia]) += C([ib,ia],[jb,ja]) * Y([jb,ja])
    * ---------------------------------------------
    */
   
   int ia = 0;
   int ja = 0;
   for(ia=0; ia < nrow_A; ia++) {
   for(ja=0; ja < ncol_A; ja++) {
         double aij = a_(ia,ja);

         int ib = 0;
         for(ib=0; ib < nrow_B; ib++) {
             int istartb = b.getRowPtr(ib);
             int iendb = b.getRowPtr(ib + 1);

             int kb = 0;
             for(kb=istartb; kb < iendb; kb++) {
                 int jb = b.getCol(kb);
                 double bij = b.getValue(kb);
                 double cij = aij * bij;

                 int ix = (isTransB) ? jb : ib;
                 int jx = (isTransA) ? ja : ia;
                 int iy = (isTransB) ? ib : jb;
                 int jy = (isTransA) ? ia : ja;

                 xout(ix,jx) +=  (cij * yin(iy,jy));
                 };
             };
         };
      };
                 
 };
   
}




void den_csr_kron_mult(
                    const char transA,
                    const char transB,
                    const PsimagLite::Matrix<double>& a_,

                    const PsimagLite::CrsMatrix<double>& b,

                    const double* yin,
                          double* xout)

{
/*
 *   -------------------------------------------------------------
 *   A in dense matrix format
 *   B in compressed sparse ROW format
 *
 *   X += kron( A, B) * Y 
 *   that can be computed as either
 *   imethod == 1
 *
 *   X(ib,ia) += (B(ib,jb) * Y(jb,ja) ) * transpose(A(ia,ja)   or
 *               BY(ib,ja) = B(ib,jb)*Y(jb,ja)
 *               BY is nrow_B by ncol_A, need   2*nnz(B)*ncolA flops               
 *
 *   X(ib,ia) +=   BY(ib,ja) * transpose(A(ia,ja)) need 2*nnz(A)*nrowB flops
 *
 *   imethod == 2
 *
 *   X(ib,ia) += B(ib,jb) * (Y(jb,ja) * transpose(A))    or
 *                YAt(jb,ia) = Y(jb,ja) * transpose(A(ia,ja))    
 *                YAt is ncolB by nrowA, need 2*nnz(A) * ncolB flops
 *
 *   X(ib,ia) += B(ib,jb) * YAt(jb,ia)  need nnz(B) * nrowA flops
 *
 *   imethod == 3
 *
 *   X += kron(A,B) * Y   by visiting all non-zero entries in A, B        
 *
 *   this is feasible only if A and B are very sparse, need nnz(A)*nnz(B) flops
 *   -------------------------------------------------------------
 */

	const int nrow_A = a_.n_row();
	const int ncol_A = a_.n_col();
	const int nrow_B = b.row();
    const int ncol_B = b.col();

 int nnz_A = den_nnz( nrow_A, ncol_A, a_);
 int nnz_B = csr_nnz(b);
 int has_work = (nnz_A >= 1) && (nnz_B >= 1);
 if (!has_work) {
     return;
     };

 double kron_nnz = 0;
 double kron_flops = 0;
 int imethod = 1;

     const int isTransA = (transA == 'T') || (transA == 't');
     const int isTransB = (transB == 'T') || (transB == 't');
    
     const int nrow_1 = (isTransA) ? ncol_A : nrow_A;
     const int ncol_2 = (isTransB) ? nrow_B : ncol_B;


 estimate_kron_cost( nrow_1,ncol_2,nnz_A, nrow_1,ncol_2,nnz_B,
                     &kron_nnz, &kron_flops, &imethod );


 den_csr_kron_mult_method( 
                    imethod,
                    transA, 
                    transB,

                    a_,

                    b,

                    yin, 
                    xout );
}

#undef A
#undef BY
#undef YAt
#undef X
#undef Y
