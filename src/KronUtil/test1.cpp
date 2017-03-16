#include "util.h"
#include "KronUtil.h"

int main()
{

    typedef PsimagLite::Vector<double>::Type VectorType;
    typedef PsimagLite::Vector<int>::Type VectorIntType;

  int nerrors = 0;
  double threshold = 0;
  int nrow_A = 0;
  int ncol_A = 0;
  int nrow_B = 0;
  int ncol_B = 0;
  int itransA = 0;
  int itransB = 0;

  for(threshold=0.01; threshold <= 1.1; threshold += 0.1) {
  for(ncol_A=1; ncol_A <= 10; ncol_A += 3 ) {
  for(nrow_A=1; nrow_A <= 10; nrow_A += 3 ) {
  for(ncol_B=1; ncol_B <= 10; ncol_B += 3 ) {
  for(nrow_B=1; nrow_B <= 10; nrow_B += 3 ) {
  for(itransA=0; itransA <= 1; itransA++) {
  for(itransB=0; itransB <= 1; itransB++) {
     char transA = (itransA == 1) ? 'T' : 'N';
     char transB = (itransB == 1) ? 'T' : 'N';


     int imethod = 0;

     int isTransA = (itransA != 0);
     int isTransB = (itransB != 0);

     int nrow_1 = (isTransA) ? ncol_A : nrow_A;
     int ncol_1 = (isTransA) ? nrow_A : ncol_A;
      
     int nrow_2 = (isTransB) ? ncol_B : nrow_B;
     int ncol_2 = (isTransB) ? nrow_B : ncol_B;
     

     int nrow_X = nrow_2;
     int ncol_X = nrow_1;

     int nrow_Y = ncol_2;
     int ncol_Y = ncol_1;

     PsimagLite::Matrix<double> a_(nrow_A, ncol_A);
     PsimagLite::Matrix<double> b_(nrow_B, ncol_B);
     PsimagLite::Matrix<double> y_(nrow_Y, ncol_Y);
	PsimagLite::MatrixNonOwned<const double> yRef(y_);

     PsimagLite::Matrix<double> x1_(nrow_X, ncol_X);
     PsimagLite::Matrix<double> x2_(nrow_X, ncol_X);
     PsimagLite::Matrix<double> x3_(nrow_X, ncol_X);

     PsimagLite::Matrix<double> sx1_(nrow_X, ncol_X);
	  PsimagLite::MatrixNonOwned<double> sx1Ref(sx1_);
     PsimagLite::Matrix<double> sx2_(nrow_X, ncol_X);
	 PsimagLite::MatrixNonOwned<double> sx2Ref(sx2_);
     PsimagLite::Matrix<double> sx3_(nrow_X, ncol_X);
	 PsimagLite::MatrixNonOwned<double> sx3Ref(sx3_);

     den_gen_matrix( nrow_A, ncol_A,  threshold, a_);
     den_gen_matrix( nrow_B, ncol_B,  threshold, b_);
     den_gen_matrix( nrow_Y, ncol_Y,  1.0, y_);

     den_zeros( nrow_X, ncol_X, x1_);
     den_zeros( nrow_X, ncol_X, x2_);
     den_zeros( nrow_X, ncol_X, x3_);


     
     den_zeros( nrow_X, ncol_X, sx1_);
     den_zeros( nrow_X, ncol_X, sx2_);
     den_zeros( nrow_X, ncol_X, sx3_);


     imethod =1;
     den_kron_mult_method( imethod,
                           transA, transB,
                           nrow_A, ncol_A, a_,
                           nrow_B, ncol_B, b_,
                           y_,
                           x1_);

     imethod = 2;
     den_kron_mult_method( imethod,
                           transA, transB,
                           nrow_A, ncol_A, a_,
                           nrow_B, ncol_B, b_,
                           y_,
                           x2_  );
     imethod = 3;
     den_kron_mult_method( imethod,
                           transA, transB,
                           nrow_A, ncol_A, a_,
                           nrow_B, ncol_B, b_,
                           y_,
                           x3_ );

     int ix = 0;
     int jx = 0;

     for(jx=0; jx < ncol_X; jx++) {
     for(ix=0; ix < nrow_X; ix++) {
        double diff12 = ABS( x1_(ix,jx) - x2_(ix,jx) );
        double diff23 = ABS( x2_(ix,jx) - x3_(ix,jx) );
        double diff31 = ABS( x3_(ix,jx) - x1_(ix,jx) );
        double diffmax = MAX( diff12, MAX( diff23, diff31) );
        const double tol = 1.0/(1000.0*1000.0*1000.0);

        int isok = (diffmax <= tol);
        if (!isok) {
           nerrors += 1;
           printf("den: itransA %d itransB %d nrow_A %d ncol_A %d nrow_B %d ncol_B %d \n",
                   itransA, itransB,    nrow_A,ncol_A,   nrow_B, ncol_B );
           printf("ix %d, jx %d, diff12 %f, diff23 %f, diff31 %f \n",
                   ix,jx,  diff12, diff23, diff31 );
           };
        };
        };


    /*
     * ------------------
     * test sparse matrix
     * ------------------
     */
     const int max_nnz_A = nrow_A * ncol_A;
     VectorType aval(max_nnz_A);
     VectorIntType acol(max_nnz_A);
     VectorIntType arowptr(nrow_A+1);

     const int max_nnz_B = nrow_B * ncol_B;
     VectorType bval(max_nnz_B);
     VectorIntType bcol(max_nnz_B);
     VectorIntType browptr(nrow_B+1);

     den2csr( nrow_A, ncol_A, a_,
              max_nnz_A,
              arowptr, acol, aval );

     den2csr( nrow_B, ncol_B, b_,
              max_nnz_B,
              browptr, bcol, bval );




     imethod =1;
     csr_kron_mult_method( 
                     imethod,
                     transA, transB,
                     nrow_A,
                     ncol_A, 
                     arowptr, 
                     acol, 
                     aval,

                     nrow_B,
                     ncol_B, 
                     browptr, 
                     bcol, 
                     bval,

                     yRef,
                     sx1Ref);


     imethod =2;
     csr_kron_mult_method( 
                     imethod,
                     transA, transB,
                     nrow_A,
                     ncol_A, 
                     arowptr, 
                     acol, 
                     aval,

                     nrow_B,
                     ncol_B, 
                     browptr, 
                     bcol, 
                     bval,

                     yRef,
                     sx2Ref);



     imethod =3;
     csr_kron_mult_method( 
                     imethod,
                     transA, transB,
                     nrow_A,
                     ncol_A, 
                     arowptr, 
                     acol, 
                     aval,

                     nrow_B,
                     ncol_B, 
                     browptr, 
                     bcol, 
                     bval,

                     yRef,
                     sx3Ref);

     for(jx=0; jx < ncol_X; jx++) {
     for(ix=0; ix < nrow_X; ix++) {
        double diff1 = ABS( x1_(ix,jx) - sx1_(ix,jx));
        double diff2 = ABS( x2_(ix,jx) - sx2_(ix,jx));
        double diff3 = ABS( x3_(ix,jx) - sx3_(ix,jx));
        double diffmax = MAX( diff1, MAX( diff2, diff3) );
        const double tol = 1.0/(1000.0*1000.0*1000.0);
        int isok = (diffmax <= tol );
        if (!isok) {
           nerrors += 1;
           printf("csr: itransA %d itransB %d nrow_A %d ncol_A %d nrow_B %d ncol_B %d \n",
                   itransA, itransB,    nrow_A,ncol_A,   nrow_B, ncol_B );
           printf("ix %d, jx %d, diff1 %f, diff2 %f, diff3 %f \n",
                   ix,jx,  diff1, diff2, diff3 );
           };
        };
        };
    

    /*
     * ---------------------
     * test generic interface
     * ---------------------
     */

     den_zeros(nrow_X,ncol_X, x1_ );
     den_zeros(nrow_X,ncol_X, sx1_ );
 
     den_kron_mult_method( imethod,
                           transA, transB,
                           nrow_A, ncol_A, a_,
                           nrow_B, ncol_B, b_,
                           y_,
                           x1_  );
     csr_kron_mult( 
                     transA, transB,
                     nrow_A,
                     ncol_A, 
                     arowptr, 
                     acol, 
                     aval,

                     nrow_B,
                     ncol_B, 
                     browptr, 
                     bcol, 
                     bval,

                     &(y_(0,0)),
                     &sx1_(0,0));
     
     for(jx=0; jx < ncol_X; jx++) {
     for(ix=0; ix < nrow_X; ix++) {
       double diff = ABS( x1_(ix,jx) - sx1_(ix,jx) );
       const double tol = 1.0/(1000.0 * 1000.0 * 1000.0);

       int isok  = (diff <= tol);
       if (!isok) {
           nerrors += 1;
           printf("nrow_A %d ncol_A %d nrow_B %d ncol_B %d \n",
                   nrow_A,ncol_A,   nrow_B, ncol_B );
           printf("ix %d, jx %d, diff %f \n", ix,jx,diff );
           };
       };
       };

    /*
     * -----------------------
     * test mixed matrix types dense and CSR
     * -----------------------
     */
     den_zeros( nrow_X, ncol_X, sx1_ );
     den_zeros( nrow_X, ncol_X, sx2_ );
     den_zeros( nrow_X, ncol_X, sx3_ );

     imethod =1;
     den_csr_kron_mult_method( 
                     imethod,
                     transA, transB,

                     nrow_A, ncol_A, a_,

                     nrow_B, ncol_B, browptr, bcol, bval,

                     y_,
                     sx1_ );


     imethod =2;
     den_csr_kron_mult_method( 
                     imethod,
                     transA, transB,

                     nrow_A, ncol_A, a_,

                     nrow_B, ncol_B, browptr, bcol, bval,

                     y_,
                     sx2_ );



     imethod =3;
     den_csr_kron_mult_method( 
                     imethod,
                     transA, transB,

                     nrow_A, ncol_A, a_,

                     nrow_B, ncol_B, browptr, bcol, bval,

                     y_,
                     sx3_ );

     for(jx=0; jx < ncol_X; jx++) {
     for(ix=0; ix < nrow_X; ix++) {
        double diff1 = ABS( x1_(ix,jx) - sx1_(ix,jx));
        double diff2 = ABS( x2_(ix,jx) - sx2_(ix,jx));
        double diff3 = ABS( x3_(ix,jx) - sx3_(ix,jx));
        double diffmax = MAX( diff1, MAX( diff2, diff3) );
        const double tol = 1.0/(1000.0*1000.0*1000.0);
        int isok = (diffmax <= tol );
        if (!isok) {
           nerrors += 1;
           printf("den_csr: itransA %d itransB %d nrow_A %d ncol_A %d nrow_B %d ncol_B %d \n",
                   itransA, itransB,    nrow_A,ncol_A,   nrow_B, ncol_B );
           printf("ix %d, jx %d, diff1 %f, diff2 %f, diff3 %f \n",
                   ix,jx,  diff1, diff2, diff3 );
           };
        };
        };





    /*
     * -----------------------
     * test mixed matrix types CSR and dense
     * -----------------------
     */
     den_zeros( nrow_X, ncol_X, sx1_ );
     den_zeros( nrow_X, ncol_X, sx2_ );
     den_zeros( nrow_X, ncol_X, sx3_ );

     imethod =1;
     csr_den_kron_mult_method( 
                     imethod,
                     transA, transB,

                     nrow_A, ncol_A, arowptr,acol,aval,
           

                     nrow_B, ncol_B, b_,

                     y_,
                     sx1_ );


     imethod =2;
     csr_den_kron_mult_method( 
                     imethod,
                     transA, transB,

                     nrow_A, ncol_A, arowptr,acol,aval,

                     nrow_B, ncol_B, b_,

                     y_,
                     sx2_ );



     imethod =3;
     csr_den_kron_mult_method( 
                     imethod,
                     transA, transB,

                     nrow_A, ncol_A, arowptr,acol,aval,

                     nrow_B, ncol_B, b_,

                     y_,
                     sx3_ );

     for(jx=0; jx < ncol_X; jx++) {
     for(ix=0; ix < nrow_X; ix++) {
        double diff1 = ABS( x1_(ix,jx) - sx1_(ix,jx));
        double diff2 = ABS( x2_(ix,jx) - sx2_(ix,jx));
        double diff3 = ABS( x3_(ix,jx) - sx3_(ix,jx));
        double diffmax = MAX( diff1, MAX( diff2, diff3) );
        const double tol = 1.0/(1000.0*1000.0*1000.0);
        int isok = (diffmax <= tol );
        if (!isok) {
           nerrors += 1;
           printf("den_csr: itransA %d itransB %d nrow_A %d ncol_A %d nrow_B %d ncol_B %d \n",
                   itransA, itransB,    nrow_A,ncol_A,   nrow_B, ncol_B );
           printf("ix %d, jx %d, diff1 %f, diff2 %f, diff3 %f \n",
                   ix,jx,  diff1, diff2, diff3 );
           };
        };
        };



   };
   };
   };
   };
   };
   };  
   };

 if (nerrors == 0) {
    printf("pass all tests\n");
    };
 return(0);
}

#undef X1
#undef X2
#undef X3
#undef SX1
#undef SX2
#undef SX3
#undef A
#undef B
#undef Y
