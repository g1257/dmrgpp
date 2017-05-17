#include "util.h"
#include "KronUtil.h"

int main()
{
  const int idebug = 0;
  int nerrors = 0;
  double thresholdA = 0;
  double thresholdB = 0;
  int nrow_A = 0;
  int ncol_A = 0;
  int nrow_B = 0;
  int ncol_B = 0;
  int itransA = 0;
  int itransB = 0;

  for(thresholdB=0; thresholdB <= 1.1; thresholdB += 0.1) {
  for(thresholdA=0; thresholdA <= 1.1; thresholdA += 0.1) {
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
	 PsimagLite::MatrixNonOwned<double> x1Ref(x1_);
     PsimagLite::Matrix<double> x2_(nrow_X, ncol_X);
	 PsimagLite::MatrixNonOwned<double> x2Ref(x2_);
     PsimagLite::Matrix<double> x3_(nrow_X, ncol_X);
	PsimagLite::MatrixNonOwned<double> x3Ref(x3_);

     PsimagLite::Matrix<double> sx1_(nrow_X, ncol_X);
	  PsimagLite::MatrixNonOwned<double> sx1Ref(sx1_);
     PsimagLite::Matrix<double> sx2_(nrow_X, ncol_X);
	 PsimagLite::MatrixNonOwned<double> sx2Ref(sx2_);
     PsimagLite::Matrix<double> sx3_(nrow_X, ncol_X);
	 PsimagLite::MatrixNonOwned<double> sx3Ref(sx3_);

     if (thresholdA == 0) {
      if (nrow_A == ncol_A) {
          /*
           * ------------------------------------
           * special case to test identity matrix
           * ------------------------------------
           */
          if (idebug >= 1) {
            printf("nrow_A=%d, identity \n",nrow_A);
            };
          den_eye( nrow_A, ncol_A, a_ );
          assert( den_is_eye(a_) );
   
          PsimagLite::CrsMatrix<double> a(a_);
          assert( csr_is_eye(a) );
          }
       else {
        if (idebug >= 1) {
           printf("nrow_A=%d,ncol_A=%d, zeros\n",nrow_A,ncol_A);
           };

           den_zeros(nrow_A,ncol_A,a_);
           assert( den_is_zeros(a_) );
           
           PsimagLite::CrsMatrix<double> a(a_);
           assert( csr_is_zeros(a) );
          };
      }
     else {
       den_gen_matrix( nrow_A, ncol_A,  thresholdA, a_);

       PsimagLite::CrsMatrix<double> a(a_);
       assert( den_is_eye(a_) == csr_is_eye(a) );
       assert( den_is_zeros(a_) == csr_is_zeros(a) );

       };


     if ((thresholdB == 0) && (nrow_B == ncol_B)) {
       /*
        * ------------------------------------
        * special case to test identity matrix
        * ------------------------------------
        */
       if (idebug >= 1) {
          printf("nrow_B=%d, identity \n",nrow_B);
          };
       den_eye( nrow_B, ncol_B, b_ );
       assert( den_is_eye(b_));

       PsimagLite::CrsMatrix<double> b(b_);
       assert( csr_is_eye(b) );
       }
     else {
       den_gen_matrix( nrow_B, ncol_B,  thresholdB, b_);

       PsimagLite::CrsMatrix<double> b(b_);
       assert( den_is_eye(b_) == csr_is_eye(b) );
       assert( den_is_zeros(b_) == csr_is_zeros(b) );
       };

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
                           a_,
                           b_,
                           yRef.getVector(),
	                       0,
                           x1Ref.getVector(),
	                       0);

     imethod = 2;
     den_kron_mult_method( imethod,
                           transA, transB,
                           a_,
                           b_,
	                       yRef.getVector(),
	                       0,
                           x2Ref.getVector(),
	                       0);
     imethod = 3;
     den_kron_mult_method( imethod,
                           transA, transB,
                           a_,
                           b_,
	                       yRef.getVector(),
	                       0,
                           x3Ref.getVector(),
	                       0);

     int ix = 0;
     int jx = 0;

     for(jx=0; jx < ncol_X; jx++) {
     for(ix=0; ix < nrow_X; ix++) {
        double diff12 = std::abs( x1_(ix,jx) - x2_(ix,jx) );
        double diff23 = std::abs( x2_(ix,jx) - x3_(ix,jx) );
        double diff31 = std::abs( x3_(ix,jx) - x1_(ix,jx) );
        double diffmax = std::max( diff12, std::max( diff23, diff31) );
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
     PsimagLite::CrsMatrix<double> a(a_);
     assert( den_is_eye(a_) == csr_is_eye(a) );
     assert( den_is_zeros(a_) == csr_is_zeros(a) );

     PsimagLite::CrsMatrix<double> b(b_);
     assert( den_is_eye(b_) == csr_is_eye(b) );
     assert( den_is_zeros(b_) == csr_is_zeros(b) );

     imethod =1;
     csr_kron_mult_method( 
                     imethod,
                     transA,
	             transB,
                     a,

                     b,

                     yRef,
                     sx1Ref);


     imethod =2;
     csr_kron_mult_method( 
                     imethod,
                     transA,
	             transB,
                     a,

                     b,

                     yRef,
                     sx2Ref);



     imethod =3;
     csr_kron_mult_method( 
                     imethod,
                     transA,
	             transB,
                     a,

                     b,

                     yRef,
                     sx3Ref);

     for(jx=0; jx < ncol_X; jx++) {
     for(ix=0; ix < nrow_X; ix++) {
        double diff1 = std::abs( x1_(ix,jx) - sx1_(ix,jx));
        double diff2 = std::abs( x2_(ix,jx) - sx2_(ix,jx));
        double diff3 = std::abs( x3_(ix,jx) - sx3_(ix,jx));
        double diffmax = std::max( diff1, std::max( diff2, diff3) );
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
 
     den_kron_mult_method(imethod,
                          transA, transB,
                           a_,
                           b_,
                           yRef.getVector(),
	                       0,
                           x1Ref.getVector(),
	                       0);
     csr_kron_mult( 
                     transA, transB,
                     a,

                     b,

	             yRef.getVector(),
                 0,
                 sx1Ref.getVector(),
                 0);
     
     for(jx=0; jx < ncol_X; jx++) {
     for(ix=0; ix < nrow_X; ix++) {
       double diff = std::abs( x1_(ix,jx) - sx1_(ix,jx) );
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
     den_csr_kron_mult( 
                     
                     transA, transB,

                      a_,

                    b,

	             yRef.getVector(),
                 0,
                 sx1Ref.getVector(),
                 0);



     for(jx=0; jx < ncol_X; jx++) {
     for(ix=0; ix < nrow_X; ix++) {
        double diff1 = std::abs( x1_(ix,jx) - sx1_(ix,jx));
        double diff2 = 0;
        double diff3 = 0;
        double diffmax = std::max( diff1, std::max( diff2, diff3) );
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


     den_zeros( nrow_X, ncol_X, sx1_ );
     den_zeros( nrow_X, ncol_X, sx2_ );
     den_zeros( nrow_X, ncol_X, sx3_ );

     imethod =1;
     den_csr_kron_mult_method( 
                     imethod,
                     transA, transB,

                      a_,

                    b,

	             yRef.getVector(),
                 0,
                 sx1Ref.getVector(),
                 0);


     imethod =2;
     den_csr_kron_mult_method( 
                     imethod,
                     transA, transB,

                      a_,

                     b,

	             yRef.getVector(),
                 0,
                 sx2Ref.getVector(),
                 0);



     imethod =3;
     den_csr_kron_mult_method( 
                     imethod,
                     transA, transB,

                     a_,

                     b,

	             yRef.getVector(),
                 0,
                 sx3Ref.getVector(),
                 0);

     for(jx=0; jx < ncol_X; jx++) {
     for(ix=0; ix < nrow_X; ix++) {
        double diff1 = std::abs( x1_(ix,jx) - sx1_(ix,jx));
        double diff2 = std::abs( x2_(ix,jx) - sx2_(ix,jx));
        double diff3 = std::abs( x3_(ix,jx) - sx3_(ix,jx));
        double diffmax = std::max( diff1, std::max( diff2, diff3) );
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

     csr_den_kron_mult( 
                     transA, transB,

                     a,
           

                     b_,

	             yRef.getVector(),
                 0,
                 sx1Ref.getVector(),
                 0);


     for(jx=0; jx < ncol_X; jx++) {
     for(ix=0; ix < nrow_X; ix++) {
        double diff1 = std::abs( x1_(ix,jx) - sx1_(ix,jx));
        double diff2 = 0;
        double diff3 = 0;
        double diffmax = std::max( diff1, std::max( diff2, diff3) );
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



     den_zeros( nrow_X, ncol_X, sx1_ );
     den_zeros( nrow_X, ncol_X, sx2_ );
     den_zeros( nrow_X, ncol_X, sx3_ );

     imethod =1;
     csr_den_kron_mult_method( 
                     imethod,
                     transA, transB,

                     a,
           

                     b_,

	             yRef.getVector(),
                 0,
                 sx1Ref.getVector(),
                 0);

     imethod =2;
     csr_den_kron_mult_method( 
                     imethod,
                     transA, transB,

                     a,

                     b_,

	             yRef.getVector(),
                 0,
                 sx2Ref.getVector(),
                 0);



     imethod =3;
     csr_den_kron_mult_method( 
                     imethod,
                     transA, transB,

                     a,

                     b_,

	             yRef.getVector(),
                 0,
                 sx3Ref.getVector(),
                 0);

     for(jx=0; jx < ncol_X; jx++) {
     for(ix=0; ix < nrow_X; ix++) {
        double diff1 = std::abs( x1_(ix,jx) - sx1_(ix,jx));
        double diff2 = std::abs( x2_(ix,jx) - sx2_(ix,jx));
        double diff3 = std::abs( x3_(ix,jx) - sx3_(ix,jx));
        double diffmax = std::max( diff1, std::max( diff2, diff3) );
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
