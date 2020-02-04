#include "util.h"
#include "KronUtil.h"

#ifndef USE_FLOAT
typedef double RealType;
#else
typedef float RealType;
#endif

int main()
{
	const RealType denseFlopDiscount = 0.2;
	const int idebug = 0;
	int nerrors = 0;
	RealType thresholdA = 0;
	RealType thresholdB = 0;
	int nrow_A = 0;
	int ncol_A = 0;
	int nrow_B = 0;
	int ncol_B = 0;
	int itransA = 0;
	int itransB = 0;

	static const bool needsPrinting = false;
	const SizeType gemmRnb = 100;
	const SizeType threadsForGemmR = 1;
	PsimagLite::GemmR<RealType> gemmR(needsPrinting, gemmRnb, threadsForGemmR);

	for(thresholdB=0; thresholdB <= 1.1; thresholdB += 0.1) {
		for(thresholdA=0; thresholdA <= 1.1; thresholdA += 0.1) {
			for(ncol_A=1; ncol_A <= 10; ncol_A += 3 ) {
				for(nrow_A=1; nrow_A <= 10; nrow_A += 3 ) {
					for(ncol_B=1; ncol_B <= 10; ncol_B += 3 ) {
						for(nrow_B=1; nrow_B <= 10; nrow_B += 3 ) {
							for(itransA=0; itransA <= 2; itransA++) {
								for(itransB=0; itransB <= 2; itransB++) {
									char transA = (itransA == 1) ? 'T' : ((itransA == 2) ? 'C' : 'N');
									char transB = (itransB == 1) ? 'T' : ((itransB == 2) ? 'C' : 'N');


									int imethod = 0;

									int isTransA = (transA == 'T');
									int isTransB = (transB == 'T');
									int isConjTransA = (transA == 'C');
									int isConjTransB = (transB == 'C');

									int nrow_1 = (isTransA || isConjTransA) ? ncol_A : nrow_A;
									int ncol_1 = (isTransA || isConjTransA) ? nrow_A : ncol_A;

									int nrow_2 = (isTransB || isConjTransB) ? ncol_B : nrow_B;
									int ncol_2 = (isTransB || isConjTransB) ? nrow_B : ncol_B;

									int nrow_C = nrow_1 * nrow_2;
									int ncol_C = ncol_1 * ncol_2;

									int nrow_X = nrow_2;
									int ncol_X = nrow_1;

									int nrow_Y = ncol_2;
									int ncol_Y = ncol_1;



									PsimagLite::Matrix<RealType> a_(nrow_A, ncol_A);
									PsimagLite::Matrix<RealType> b_(nrow_B, ncol_B);

									PsimagLite::Matrix<RealType> y_(nrow_Y, ncol_Y);
									PsimagLite::MatrixNonOwned<const RealType> yRef(y_);

									PsimagLite::Matrix<RealType> x1_(nrow_X, ncol_X);
									PsimagLite::MatrixNonOwned<RealType> x1Ref(x1_);
									PsimagLite::Matrix<RealType> x2_(nrow_X, ncol_X);
									PsimagLite::MatrixNonOwned<RealType> x2Ref(x2_);
									PsimagLite::Matrix<RealType> x3_(nrow_X, ncol_X);
									PsimagLite::MatrixNonOwned<RealType> x3Ref(x3_);
									PsimagLite::Matrix<RealType> x4_(nrow_X, ncol_X);
									PsimagLite::MatrixNonOwned<RealType> x4Ref(x4_);

									PsimagLite::Matrix<RealType> sx1_(nrow_X, ncol_X);
									PsimagLite::MatrixNonOwned<RealType> sx1Ref(sx1_);
									PsimagLite::Matrix<RealType> sx2_(nrow_X, ncol_X);
									PsimagLite::MatrixNonOwned<RealType> sx2Ref(sx2_);
									PsimagLite::Matrix<RealType> sx3_(nrow_X, ncol_X);
									PsimagLite::MatrixNonOwned<RealType> sx3Ref(sx3_);
									PsimagLite::Matrix<RealType> sx4_(nrow_X, ncol_X);
									PsimagLite::MatrixNonOwned<RealType> sx4Ref(sx4_);

									if (thresholdA == 0) {
										if (nrow_A == ncol_A) {
											/*
		   * ------------------------------------
		   * special case to test identity matrix
		   * ------------------------------------
		   */
											if (idebug >= 1) {
												printf("nrow_A=%d, identity \n",nrow_A);
											}
											den_eye( nrow_A, ncol_A, a_ );
											assert( den_is_eye(a_) );

											PsimagLite::CrsMatrix<RealType> a(a_);
											assert( csr_is_eye(a) );
										}
										else {
											if (idebug >= 1) {
												printf("nrow_A=%d,ncol_A=%d, zeros\n",nrow_A,ncol_A);
											}

											den_zeros(nrow_A,ncol_A,a_);
											assert( den_is_zeros(a_) );

											PsimagLite::CrsMatrix<RealType> a(a_);
											assert( csr_is_zeros(a) );
										}
									}
									else {
										den_gen_matrix( nrow_A, ncol_A,  thresholdA, a_);

										PsimagLite::CrsMatrix<RealType> a(a_);
										assert( den_is_eye(a_) == csr_is_eye(a) );
										assert( den_is_zeros(a_) == csr_is_zeros(a) );

									}


									if ((thresholdB == 0) && (nrow_B == ncol_B)) {
										/*
		* ------------------------------------
		* special case to test identity matrix
		* ------------------------------------
		*/
										if (idebug >= 1) {
											printf("nrow_B=%d, identity \n",nrow_B);
										}
										den_eye( nrow_B, ncol_B, b_ );
										assert( den_is_eye(b_));

										PsimagLite::CrsMatrix<RealType> b(b_);
										assert( csr_is_eye(b) );
									}
									else {
										den_gen_matrix( nrow_B, ncol_B,  thresholdB, b_);

										PsimagLite::CrsMatrix<RealType> b(b_);
										assert( den_is_eye(b_) == csr_is_eye(b) );
										assert( den_is_zeros(b_) == csr_is_zeros(b) );
									}

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
									                      0,
									                      gemmR);

									imethod = 2;
									den_kron_mult_method( imethod,
									                      transA, transB,
									                      a_,
									                      b_,
									                      yRef.getVector(),
									                      0,
									                      x2Ref.getVector(),
									                      0,
									                      gemmR);
									imethod = 3;
									den_kron_mult_method( imethod,
									                      transA, transB,
									                      a_,
									                      b_,
									                      yRef.getVector(),
									                      0,
									                      x3Ref.getVector(),
									                      0,
									                      gemmR);

									// ------------------
									// form C = kron(A,B)
									// ------------------
									PsimagLite::Matrix<RealType> c_(nrow_C, ncol_C);

									den_kron_form_general(
									            transA, transB,
									            nrow_A, ncol_A, a_,
									            nrow_B, ncol_B, b_,
									            c_ );

									// -----------------------
									// perform matrix-multiply
									// -----------------------
									{
										const char trans1 = 'N';
										const char trans2 = 'N';
										const RealType alpha = 1.0;
										const RealType beta = 0.0;

										// ------------------------------
										// reshape X, Y as column vectors
										// ------------------------------
										const int mm = nrow_X*ncol_X;
										const int nn = 1;
										const int kk = ncol_C;


										const int ld1 = nrow_C;
										const int ld2 = nrow_Y*ncol_Y;
										const int ld3 = nrow_X*ncol_X;

										const RealType * const pA = &(c_(0,0));
										const RealType * const pB = &(yRef.getVector()[0]);
										RealType *pC = &(x4Ref.getVector()[0]);
										psimag::BLAS::GEMM( trans1, trans2, mm,nn,kk,
										                    alpha, pA, ld1,  pB, ld2,
										                    beta,  pC, ld3 );
									}




									int ix = 0;
									int jx = 0;

									for(jx=0; jx < ncol_X; jx++) {
										for(ix=0; ix < nrow_X; ix++) {
											RealType diff12 = std::abs( x1_(ix,jx) - x2_(ix,jx) );
											RealType diff23 = std::abs( x2_(ix,jx) - x3_(ix,jx) );
											RealType diff31 = std::abs( x3_(ix,jx) - x1_(ix,jx) );
											RealType diff41 = std::abs( x4_(ix,jx) - x1_(ix,jx) );
											RealType diffmax = std::max( diff41,
											                             std::max( diff12, std::max( diff23, diff31) ) );
											const RealType tol = 1.0/(1000.0*1000.0*1000.0);

											int isok = (diffmax <= tol);
											if (!isok) {
												nerrors += 1;
												printf("den: transA=%c, itransA %d transB=%c, itransB %d nrow_A %d ncol_A %d nrow_B %d ncol_B %d \n",
												       transA,itransA, transB,itransB,    nrow_A,ncol_A,   nrow_B, ncol_B );
												printf("ix %d, jx %d, diff12 %f, diff23 %f, diff31 %f diff41 %f \n",
												       ix,    jx,    diff12,    diff23,    diff31,   diff41 );
											}
										}
									}


									/*
	 * ------------------
	 * test sparse matrix
	 * ------------------
	 */
									PsimagLite::CrsMatrix<RealType> a(a_);
									assert( den_is_eye(a_) == csr_is_eye(a) );
									assert( den_is_zeros(a_) == csr_is_zeros(a) );

									PsimagLite::CrsMatrix<RealType> b(b_);
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
											RealType diff1 = std::abs( x1_(ix,jx) - sx1_(ix,jx));
											RealType diff2 = std::abs( x2_(ix,jx) - sx2_(ix,jx));
											RealType diff3 = std::abs( x3_(ix,jx) - sx3_(ix,jx));
											RealType diffmax = std::max( diff1, std::max( diff2, diff3) );
											const RealType tol = 1.0/(1000.0*1000.0*1000.0);
											int isok = (diffmax <= tol );
											if (!isok) {
												nerrors += 1;
												printf("csr: itransA %d itransB %d nrow_A %d ncol_A %d nrow_B %d ncol_B %d \n",
												       itransA, itransB,    nrow_A,ncol_A,   nrow_B, ncol_B );
												printf("ix %d, jx %d, diff1 %f, diff2 %f, diff3 %f \n",
												       ix,jx,  diff1, diff2, diff3 );
											}
										}
									}


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
									                     0,
									                     gemmR);

									csr_kron_mult(transA,
									              transB,
									              a,
									              b,
									              yRef.getVector(),
									              0,
									              sx1Ref.getVector(),
									              0,
									              denseFlopDiscount);

									for(jx=0; jx < ncol_X; jx++) {
										for(ix=0; ix < nrow_X; ix++) {
											RealType diff = std::abs( x1_(ix,jx) - sx1_(ix,jx) );
											const RealType tol = 1.0/(1000.0 * 1000.0 * 1000.0);

											int isok  = (diff <= tol);
											if (!isok) {
												nerrors += 1;
												printf("nrow_A %d ncol_A %d nrow_B %d ncol_B %d \n",
												       nrow_A,ncol_A,   nrow_B, ncol_B );
												printf("ix %d, jx %d, diff %f \n", ix,jx,diff );
											}
										}
									}

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
									            0,
									            denseFlopDiscount,
									            gemmR);



									for(jx=0; jx < ncol_X; jx++) {
										for(ix=0; ix < nrow_X; ix++) {
											RealType diff1 = std::abs( x1_(ix,jx) - sx1_(ix,jx));
											RealType diff2 = 0;
											RealType diff3 = 0;
											RealType diffmax = std::max( diff1, std::max( diff2, diff3) );
											const RealType tol = 1.0/(1000.0*1000.0*1000.0);
											int isok = (diffmax <= tol );
											if (!isok) {
												nerrors += 1;
												printf("den_csr: itransA %d itransB %d nrow_A %d ncol_A %d nrow_B %d ncol_B %d \n",
												       itransA, itransB,    nrow_A,ncol_A,   nrow_B, ncol_B );
												printf("ix %d, jx %d, diff1 %f, diff2 %f, diff3 %f \n",
												       ix,jx,  diff1, diff2, diff3 );
											}
										}
									}


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
									            0,
									            gemmR);


									imethod =2;
									den_csr_kron_mult_method(
									            imethod,
									            transA, transB,

									            a_,

									            b,

									            yRef.getVector(),
									            0,
									            sx2Ref.getVector(),
									            0,
									            gemmR);



									imethod =3;
									den_csr_kron_mult_method(imethod,
									                         transA,
									                         transB,
									                         a_,
									                         b,
									                         yRef.getVector(),
									                         0,
									                         sx3Ref.getVector(),
									                         0,
									                         gemmR);

									for(jx=0; jx < ncol_X; jx++) {
										for(ix=0; ix < nrow_X; ix++) {
											RealType diff1 = std::abs( x1_(ix,jx) - sx1_(ix,jx));
											RealType diff2 = std::abs( x2_(ix,jx) - sx2_(ix,jx));
											RealType diff3 = std::abs( x3_(ix,jx) - sx3_(ix,jx));
											RealType diffmax = std::max( diff1, std::max( diff2, diff3) );
											const RealType tol = 1.0/(1000.0*1000.0*1000.0);
											int isok = (diffmax <= tol );
											if (!isok) {
												nerrors += 1;
												printf("den_csr: itransA %d itransB %d nrow_A %d ncol_A %d nrow_B %d ncol_B %d \n",
												       itransA, itransB,    nrow_A,ncol_A,   nrow_B, ncol_B );
												printf("ix %d, jx %d, diff1 %f, diff2 %f, diff3 %f \n",
												       ix,jx,  diff1, diff2, diff3 );
											}
										}
									}





									/*
	 * -----------------------
	 * test mixed matrix types CSR and dense
	 * -----------------------
	 */
									den_zeros( nrow_X, ncol_X, sx1_ );

									csr_den_kron_mult(transA,
									                  transB,
									                  a,
									                  b_,
									                  yRef.getVector(),
									                  0,
									                  sx1Ref.getVector(),
									                  0,
									                  denseFlopDiscount,
									                  gemmR);


									for(jx=0; jx < ncol_X; jx++) {
										for(ix=0; ix < nrow_X; ix++) {
											RealType diff1 = std::abs( x1_(ix,jx) - sx1_(ix,jx));
											RealType diff2 = 0;
											RealType diff3 = 0;
											RealType diffmax = std::max( diff1, std::max( diff2, diff3) );
											const RealType tol = 1.0/(1000.0*1000.0*1000.0);
											int isok = (diffmax <= tol );
											if (!isok) {
												nerrors += 1;
												printf("den_csr: itransA %d itransB %d nrow_A %d ncol_A %d nrow_B %d ncol_B %d \n",
												       itransA, itransB,    nrow_A,ncol_A,   nrow_B, ncol_B );
												printf("ix %d, jx %d, diff1 %f, diff2 %f, diff3 %f \n",
												       ix,jx,  diff1, diff2, diff3 );
											}
										}
									}



									den_zeros( nrow_X, ncol_X, sx1_ );
									den_zeros( nrow_X, ncol_X, sx2_ );
									den_zeros( nrow_X, ncol_X, sx3_ );

									imethod =1;
									csr_den_kron_mult_method(imethod,
									                         transA,
									                         transB,
									                         a,
									                         b_,
									                         yRef.getVector(),
									                         0,
									                         sx1Ref.getVector(),
									                         0,
									                         gemmR);

									imethod =2;
									csr_den_kron_mult_method(imethod,
									                         transA,
									                         transB,
									                         a,
									                         b_,
									                         yRef.getVector(),
									                         0,
									                         sx2Ref.getVector(),
									                         0,
									                         gemmR);



									imethod =3;
									csr_den_kron_mult_method(imethod,
									                         transA,
									                         transB,
									                         a,
									                         b_,
									                         yRef.getVector(),
									                         0,
									                         sx3Ref.getVector(),
									                         0,
									                         gemmR);

									for(jx=0; jx < ncol_X; jx++) {
										for(ix=0; ix < nrow_X; ix++) {
											RealType diff1 = std::abs( x1_(ix,jx) - sx1_(ix,jx));
											RealType diff2 = std::abs( x2_(ix,jx) - sx2_(ix,jx));
											RealType diff3 = std::abs( x3_(ix,jx) - sx3_(ix,jx));
											RealType diffmax = std::max( diff1, std::max( diff2, diff3) );
											const RealType tol = 1.0/(1000.0*1000.0*1000.0);
											int isok = (diffmax <= tol );
											if (!isok) {
												nerrors += 1;
												printf("den_csr: itransA %d itransB %d nrow_A %d ncol_A %d nrow_B %d ncol_B %d \n",
												       itransA, itransB,    nrow_A,ncol_A,   nrow_B, ncol_B );
												printf("ix %d, jx %d, diff1 %f, diff2 %f, diff3 %f \n",
												       ix,jx,  diff1, diff2, diff3 );
											}
										}
									}



								}
							}
						}
					}
				}
			}
		}
	}

	if (nerrors == 0) {
		printf("pass all tests\n");
	}
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
