#include "KronUtil.h"
#include "util.h"

#include <Kokkos_Core.hpp>
#define CATCH_CONFIG_RUNNER
#include <catch2/catch_session.hpp>
#include <catch2/catch_test_macros.hpp>
#include <complex>

template <typename T> void run_kron_checks()
{
	using RealT                   = typename PsimagLite::Real<T>::Type;
	const RealT denseFlopDiscount = RealT(0.2);
	double      thresholdA        = 0.0;
	double      thresholdB        = 0.0;
	int         nrow_A            = 0;
	int         ncol_A            = 0;
	int         nrow_B            = 0;
	int         ncol_B            = 0;
	int         itransA           = 0;
	int         itransB           = 0;

	static const bool    needsPrinting   = false;
	const SizeType       gemmRnb         = 49;
	const SizeType       threadsForGemmR = 1;
	PsimagLite::GemmR<T> gemmR(needsPrinting, gemmRnb, threadsForGemmR);

	for (thresholdB = 0; thresholdB <= 1.1; thresholdB += 0.1) {
		for (thresholdA = 0; thresholdA <= 1.1; thresholdA += 0.1) {
			for (ncol_A = 1; ncol_A <= 7; ncol_A += 3) {
				for (nrow_A = 1; nrow_A <= 7; nrow_A += 3) {
					for (ncol_B = 1; ncol_B <= 7; ncol_B += 3) {
						for (nrow_B = 1; nrow_B <= 10; nrow_B += 3) {
							for (itransA = 0; itransA <= 2; itransA++) {
								for (itransB = 0; itransB <= 2;
								     itransB++) {
									char transA = (itransA == 1)
									    ? 'T'
									    : ((itransA == 2)
									           ? 'C'
									           : 'N');
									char transB = (itransB == 1)
									    ? 'T'
									    : ((itransB == 2)
									           ? 'C'
									           : 'N');

									int imethod = 0;

									int isTransA
									    = (transA == 'T');
									int isTransB
									    = (transB == 'T');
									int isConjTransA
									    = (transA == 'C');
									int isConjTransB
									    = (transB == 'C');

									int nrow_1
									    = (isTransA
									       || isConjTransA)
									    ? ncol_A
									    : nrow_A;
									int ncol_1
									    = (isTransA
									       || isConjTransA)
									    ? nrow_A
									    : ncol_A;

									int nrow_2
									    = (isTransB
									       || isConjTransB)
									    ? ncol_B
									    : nrow_B;
									int ncol_2
									    = (isTransB
									       || isConjTransB)
									    ? nrow_B
									    : ncol_B;

									int nrow_C
									    = nrow_1 * nrow_2;
									int ncol_C
									    = ncol_1 * ncol_2;

									int nrow_X = nrow_2;
									int ncol_X = nrow_1;

									int nrow_Y = ncol_2;
									int ncol_Y = ncol_1;

									PsimagLite::Matrix<T> a_(
									    nrow_A, ncol_A);
									PsimagLite::Matrix<T> b_(
									    nrow_B, ncol_B);

									PsimagLite::Matrix<T> y_(
									    nrow_Y, ncol_Y);
									PsimagLite::MatrixNonOwned<
									    const T>
									    yRef(y_);

									PsimagLite::Matrix<T> x1_(
									    nrow_X, ncol_X);
									PsimagLite::MatrixNonOwned<
									    T>
									    x1Ref(x1_);
									PsimagLite::Matrix<T> x2_(
									    nrow_X, ncol_X);
									PsimagLite::MatrixNonOwned<
									    T>
									    x2Ref(x2_);
									PsimagLite::Matrix<T> x3_(
									    nrow_X, ncol_X);
									PsimagLite::MatrixNonOwned<
									    T>
									    x3Ref(x3_);
									PsimagLite::Matrix<T> x4_(
									    nrow_X, ncol_X);
									PsimagLite::MatrixNonOwned<
									    T>
									    x4Ref(x4_);

									PsimagLite::Matrix<T> sx1_(
									    nrow_X, ncol_X);
									PsimagLite::MatrixNonOwned<
									    T>
									    sx1Ref(sx1_);
									PsimagLite::Matrix<T> sx2_(
									    nrow_X, ncol_X);
									PsimagLite::MatrixNonOwned<
									    T>
									    sx2Ref(sx2_);
									PsimagLite::Matrix<T> sx3_(
									    nrow_X, ncol_X);
									PsimagLite::MatrixNonOwned<
									    T>
									    sx3Ref(sx3_);
									PsimagLite::Matrix<T> sx4_(
									    nrow_X, ncol_X);
									PsimagLite::MatrixNonOwned<
									    T>
									    sx4Ref(sx4_);

									if (thresholdA == 0) {
										if (nrow_A
										    == ncol_A) {
											den_eye(
											    nrow_A,
											    ncol_A,
											    a_);
											REQUIRE(
											    den_is_eye(
											        a_));
											PsimagLite::
											    CrsMatrix<
											        T>
											        a(a_);
											REQUIRE(
											    csr_is_eye(
											        a));
										} else {
											den_zeros(
											    nrow_A,
											    ncol_A,
											    a_);
											REQUIRE(
											    den_is_zeros(
											        a_));
											PsimagLite::
											    CrsMatrix<
											        T>
											        a(a_);
											REQUIRE(
											    csr_is_zeros(
											        a));
										}
									} else {
										den_gen_matrix(
										    nrow_A,
										    ncol_A,
										    thresholdA,
										    a_);
										PsimagLite::
										    CrsMatrix<T>
										        a(a_);
										REQUIRE(
										    den_is_eye(a_)
										    == csr_is_eye(
										        a));
										REQUIRE(
										    den_is_zeros(a_)
										    == csr_is_zeros(
										        a));
									}

									if ((thresholdB == 0)
									    && (nrow_B == ncol_B)) {
										den_eye(nrow_B,
										        ncol_B,
										        b_);
										REQUIRE(
										    den_is_eye(b_));
										PsimagLite::
										    CrsMatrix<T>
										        b(b_);
										REQUIRE(
										    csr_is_eye(b));
									} else {
										den_gen_matrix(
										    nrow_B,
										    ncol_B,
										    thresholdB,
										    b_);
										PsimagLite::
										    CrsMatrix<T>
										        b(b_);
										REQUIRE(
										    den_is_eye(b_)
										    == csr_is_eye(
										        b));
										REQUIRE(
										    den_is_zeros(b_)
										    == csr_is_zeros(
										        b));
									}

									den_gen_matrix(nrow_Y,
									               ncol_Y,
									               1.0,
									               y_);

									den_zeros(
									    nrow_X, ncol_X, x1_);
									den_zeros(
									    nrow_X, ncol_X, x2_);
									den_zeros(
									    nrow_X, ncol_X, x3_);

									den_zeros(
									    nrow_X, ncol_X, sx1_);
									den_zeros(
									    nrow_X, ncol_X, sx2_);
									den_zeros(
									    nrow_X, ncol_X, sx3_);

									imethod = 1;
									den_kron_mult_method(
									    imethod,
									    transA,
									    transB,
									    a_,
									    b_,
									    yRef.getVector(),
									    0,
									    x1Ref.getVector(),
									    0,
									    gemmR);

									imethod = 2;
									den_kron_mult_method(
									    imethod,
									    transA,
									    transB,
									    a_,
									    b_,
									    yRef.getVector(),
									    0,
									    x2Ref.getVector(),
									    0,
									    gemmR);
									imethod = 3;
									den_kron_mult_method(
									    imethod,
									    transA,
									    transB,
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
									PsimagLite::Matrix<T> c_(
									    nrow_C, ncol_C);

									den_kron_form_general(
									    transA,
									    transB,
									    nrow_A,
									    ncol_A,
									    a_,
									    nrow_B,
									    ncol_B,
									    b_,
									    c_);

									// -----------------------
									// perform matrix-multiply
									// -----------------------
									{
										const char trans1
										    = 'N';
										const char trans2
										    = 'N';
										const T alpha = 1.0;
										const T beta  = 0.0;

										// ------------------------------
										// reshape X, Y as
										// column vectors
										// ------------------------------
										const int mm
										    = nrow_X
										    * ncol_X;
										const int nn = 1;
										const int kk
										    = ncol_C;

										const int ld1
										    = nrow_C;
										const int ld2
										    = nrow_Y
										    * ncol_Y;
										const int ld3
										    = nrow_X
										    * ncol_X;

										const T* const pA
										    = &(c_(0, 0));
										const T* const pB = &(
										    yRef.getVector()
										        [0]);
										T* pC = &(
										    x4Ref
										        .getVector()
										            [0]);
										psimag::BLAS::GEMM(
										    trans1,
										    trans2,
										    mm,
										    nn,
										    kk,
										    alpha,
										    pA,
										    ld1,
										    pB,
										    ld2,
										    beta,
										    pC,
										    ld3);
									}

									for (int jx = 0;
									     jx < ncol_X;
									     ++jx) {
										for (int ix = 0;
										     ix < nrow_X;
										     ++ix) {
											auto diff12 = std::abs(
											    x1_(ix,
											        jx)
											    - x2_(
											        ix,
											        jx));
											auto diff23 = std::abs(
											    x2_(ix,
											        jx)
											    - x3_(
											        ix,
											        jx));
											auto diff31 = std::abs(
											    x3_(ix,
											        jx)
											    - x1_(
											        ix,
											        jx));
											auto diff41 = std::abs(
											    x4_(ix,
											        jx)
											    - x1_(
											        ix,
											        jx));
											auto diffmax = std::max(
											    diff41,
											    std::max(
											        diff12,
											        std::max(
											            diff23,
											            diff31)));
											const double
											    tol
											    = 1.0
											    / (1000.0
											       * 1000.0
											       * 1000.0);
											if (diffmax
											    > tol) {
												INFO(
												    "den: "
												    "transA"
												    "="
												    << transA
												    << ", "
												       "itr"
												       "ans"
												       "A "
												    << itransA
												    << " tr"
												       "ans"
												       "B="
												    << transB
												    << ", "
												       "itr"
												       "ans"
												       "B "
												    << itransB
												    << " nr"
												       "ow_"
												       "A "
												    << nrow_A
												    << " nc"
												       "ol_"
												       "A "
												    << ncol_A
												    << " nr"
												       "ow_"
												       "B "
												    << nrow_B
												    << " nco"
												       "l_"
												       "B "
												    << ncol_B
												    << '\n'
												    << "ix "
												    << ix
												    << ", "
												       "jx "
												    << jx
												    << ", "
												       "dif"
												       "f12"
												       " "
												    << diff12
												    << ", "
												       "dif"
												       "f23"
												       " "
												    << diff23
												    << ", "
												       "dif"
												       "f31"
												       " "
												    << diff31
												    << " di"
												       "ff4"
												       "1 "
												    << diff41
												    << '\n');
												REQUIRE(
												    diffmax
												    <= tol);
											}
										}
									}

									/*
									 * ------------------
									 * test sparse matrix
									 * ------------------
									 */
									PsimagLite::CrsMatrix<T> a(
									    a_);
									REQUIRE(den_is_eye(a_)
									        == csr_is_eye(a));
									REQUIRE(den_is_zeros(a_)
									        == csr_is_zeros(a));
									PsimagLite::CrsMatrix<T> b(
									    b_);
									REQUIRE(den_is_eye(b_)
									        == csr_is_eye(b));
									REQUIRE(den_is_zeros(b_)
									        == csr_is_zeros(b));

									imethod = 1;
									csr_kron_mult_method(
									    imethod,
									    transA,
									    transB,
									    a,

									    b,

									    yRef,
									    sx1Ref);

									imethod = 2;
									csr_kron_mult_method(
									    imethod,
									    transA,
									    transB,
									    a,

									    b,

									    yRef,
									    sx2Ref);

									imethod = 3;
									csr_kron_mult_method(
									    imethod,
									    transA,
									    transB,
									    a,

									    b,

									    yRef,
									    sx3Ref);

									for (int jx = 0;
									     jx < ncol_X;
									     ++jx) {
										for (int ix = 0;
										     ix < nrow_X;
										     ++ix) {
											auto diff1 = std::abs(
											    x1_(ix,
											        jx)
											    - sx1_(
											        ix,
											        jx));
											auto diff2 = std::abs(
											    x2_(ix,
											        jx)
											    - sx2_(
											        ix,
											        jx));
											auto diff3 = std::abs(
											    x3_(ix,
											        jx)
											    - sx3_(
											        ix,
											        jx));
											auto diffmax = std::max(
											    diff1,
											    std::max(
											        diff2,
											        diff3));
											const double
											    tol
											    = 1.0
											    / (1000.0
											       * 1000.0
											       * 1000.0);
											if (diffmax
											    > tol) {
												INFO(
												    "csr: "
												    "transA"
												    "="
												    << transA
												    << ", "
												       "itr"
												       "ans"
												       "A "
												    << itransA
												    << " tr"
												       "ans"
												       "B="
												    << transB
												    << ", "
												       "itr"
												       "ans"
												       "B "
												    << itransB
												    << " nr"
												       "ow_"
												       "A "
												    << nrow_A
												    << " nc"
												       "ol_"
												       "A "
												    << ncol_A
												    << " nr"
												       "ow_"
												       "B "
												    << nrow_B
												    << " nco"
												       "l_"
												       "B "
												    << ncol_B
												    << '\n'
												    << "ix "
												    << ix
												    << ", "
												       "jx "
												    << jx
												    << ", "
												       "dif"
												       "f1 "
												    << diff1
												    << ", "
												       "dif"
												       "f2 "
												    << diff2
												    << ", "
												       "dif"
												       "f3 "
												    << diff3
												    << '\n');
												REQUIRE(
												    diffmax
												    <= tol);
											}
										}
									}

									/*
									 * ---------------------
									 * test generic interface
									 * ---------------------
									 */

									den_zeros(
									    nrow_X, ncol_X, x1_);
									den_zeros(
									    nrow_X, ncol_X, sx1_);

									den_kron_mult_method(
									    imethod,
									    transA,
									    transB,
									    a_,
									    b_,
									    yRef.getVector(),
									    0,
									    x1Ref.getVector(),
									    0,
									    gemmR);

									csr_kron_mult(
									    transA,
									    transB,
									    a,
									    b,
									    yRef.getVector(),
									    0,
									    sx1Ref.getVector(),
									    0,
									    RealT(
									        denseFlopDiscount));

									for (int jx = 0;
									     jx < ncol_X;
									     ++jx) {
										for (int ix = 0;
										     ix < nrow_X;
										     ++ix) {
											auto diff = std::abs(
											    x1_(ix,
											        jx)
											    - sx1_(
											        ix,
											        jx));
											const double
											    tol
											    = 1.0
											    / (1000.0
											       * 1000.0
											       * 1000.0);
											if (diff
											    > tol) {
												INFO(
												    "nrow_A "
												    << nrow_A
												    << " ncol_A "
												    << ncol_A
												    << " nrow_B "
												    << nrow_B
												    << " ncol_B "
												    << ncol_B
												    << '\n'
												    << "ix "
												    << ix
												    << ", jx "
												    << jx
												    << ", diff "
												    << diff
												    << '\n');
												REQUIRE(
												    diff
												    <= tol);
											}
										}
									}

									/*
									 * -----------------------
									 * test mixed matrix types
									 * dense and CSR
									 * -----------------------
									 */

									den_zeros(
									    nrow_X, ncol_X, sx1_);
									den_csr_kron_mult(

									    transA,
									    transB,

									    a_,

									    b,

									    yRef.getVector(),
									    0,
									    sx1Ref.getVector(),
									    0,
									    RealT(
									        denseFlopDiscount),
									    gemmR);

									for (int jx = 0;
									     jx < ncol_X;
									     ++jx) {
										for (int ix = 0;
										     ix < nrow_X;
										     ++ix) {
											auto diff1 = std::abs(
											    x1_(ix,
											        jx)
											    - sx1_(
											        ix,
											        jx));
											auto diff2
											    = 0.0;
											auto diff3
											    = 0.0;
											auto diffmax = std::max(
											    diff1,
											    std::max(
											        diff2,
											        diff3));
											const double
											    tol
											    = 1.0
											    / (1000.0
											       * 1000.0
											       * 1000.0);
											if (diffmax
											    > tol) {
												INFO(
												    "den_csr: "
												    "itr"
												    "ans"
												    "A "
												    << itransA
												    << "itr"
												       "ans"
												       "B "
												    << itransB
												    << " nr"
												       "ow_"
												       "A "
												    << nrow_A
												    << " nc"
												       "ol_"
												       "A "
												    << ncol_A
												    << " nr"
												       "ow_"
												       "B "
												    << nrow_B
												    << " nco"
												       "l_"
												       "B "
												    << ncol_B
												    << '\n'
												    << "ix "
												    << ix
												    << ", "
												       "jx "
												    << jx
												    << ", "
												       "dif"
												       "f1 "
												    << diff1
												    << ", "
												       "dif"
												       "f2 "
												    << diff2
												    << ", "
												       "dif"
												       "f3 "
												    << diff3
												    << '\n');
												REQUIRE(
												    diffmax
												    <= tol);
											}
										}
									}

									den_zeros(
									    nrow_X, ncol_X, sx1_);
									den_zeros(
									    nrow_X, ncol_X, sx2_);
									den_zeros(
									    nrow_X, ncol_X, sx3_);

									imethod = 1;
									den_csr_kron_mult_method(
									    imethod,
									    transA,
									    transB,

									    a_,

									    b,

									    yRef.getVector(),
									    0,
									    sx1Ref.getVector(),
									    0,
									    gemmR);

									imethod = 2;
									den_csr_kron_mult_method(
									    imethod,
									    transA,
									    transB,

									    a_,

									    b,

									    yRef.getVector(),
									    0,
									    sx2Ref.getVector(),
									    0,
									    gemmR);

									imethod = 3;
									den_csr_kron_mult_method(
									    imethod,
									    transA,
									    transB,
									    a_,
									    b,
									    yRef.getVector(),
									    0,
									    sx3Ref.getVector(),
									    0,
									    gemmR);

									for (int jx = 0;
									     jx < ncol_X;
									     ++jx) {
										for (int ix = 0;
										     ix < nrow_X;
										     ++ix) {
											auto diff1 = std::abs(
											    x1_(ix,
											        jx)
											    - sx1_(
											        ix,
											        jx));
											auto diff2 = std::abs(
											    x2_(ix,
											        jx)
											    - sx2_(
											        ix,
											        jx));
											auto diff3 = std::abs(
											    x3_(ix,
											        jx)
											    - sx3_(
											        ix,
											        jx));
											auto diffmax = std::max(
											    diff1,
											    std::max(
											        diff2,
											        diff3));
											const double
											    tol
											    = 1.0
											    / (1000.0
											       * 1000.0
											       * 1000.0);
											if (diffmax
											    > tol) {
												INFO(
												    "den_csr: "
												    "itr"
												    "ans"
												    "A "
												    << itransA
												    << "itr"
												       "ans"
												       "B "
												    << itransB
												    << " nr"
												       "ow_"
												       "A "
												    << nrow_A
												    << " nc"
												       "ol_"
												       "A "
												    << ncol_A
												    << " nr"
												       "ow_"
												       "B "
												    << nrow_B
												    << " nco"
												       "l_"
												       "B "
												    << ncol_B
												    << '\n'
												    << "ix "
												    << ix
												    << ", "
												       "jx "
												    << jx
												    << ", "
												       "dif"
												       "f1 "
												    << diff1
												    << ", "
												       "dif"
												       "f2 "
												    << diff2
												    << ", "
												       "dif"
												       "f3 "
												    << diff3
												    << '\n');
												REQUIRE(
												    diffmax
												    <= tol);
											}
										}
									}

									/*
									 * -----------------------
									 * test mixed matrix types
									 * CSR and dense
									 * -----------------------
									 */
									den_zeros(
									    nrow_X, ncol_X, sx1_);

									csr_den_kron_mult(
									    transA,
									    transB,
									    a,
									    b_,
									    yRef.getVector(),
									    0,
									    sx1Ref.getVector(),
									    0,
									    RealT(
									        denseFlopDiscount),
									    gemmR);

									for (int jx = 0;
									     jx < ncol_X;
									     ++jx) {
										for (int ix = 0;
										     ix < nrow_X;
										     ++ix) {
											auto diff1 = std::abs(
											    x1_(ix,
											        jx)
											    - sx1_(
											        ix,
											        jx));
											auto diff2
											    = 0.0;
											auto diff3
											    = 0.0;
											auto diffmax = std::max(
											    diff1,
											    std::max(
											        diff2,
											        diff3));
											const double
											    tol
											    = 1.0
											    / (1000.0
											       * 1000.0
											       * 1000.0);
											if (diffmax
											    > tol) {
												INFO(
												    "den_csr: "
												    "itr"
												    "ans"
												    "A "
												    << itransA
												    << "itr"
												       "ans"
												       "B "
												    << itransB
												    << " nr"
												       "ow_"
												       "A "
												    << nrow_A
												    << " nc"
												       "ol_"
												       "A "
												    << ncol_A
												    << " nr"
												       "ow_"
												       "B "
												    << nrow_B
												    << " nco"
												       "l_"
												       "B "
												    << ncol_B
												    << '\n'
												    << "ix "
												    << ix
												    << ", "
												       "jx "
												    << jx
												    << ", "
												       "dif"
												       "f1 "
												    << diff1
												    << ", "
												       "dif"
												       "f2 "
												    << diff2
												    << ", "
												       "dif"
												       "f3 "
												    << diff3
												    << '\n');
												REQUIRE(
												    diffmax
												    <= tol);
											}
										}
									}

									den_zeros(
									    nrow_X, ncol_X, sx1_);
									den_zeros(
									    nrow_X, ncol_X, sx2_);
									den_zeros(
									    nrow_X, ncol_X, sx3_);

									imethod              = 1;
									const SizeType izero = 0;
									csr_den_kron_mult_method(
									    imethod,
									    transA,
									    transB,
									    a,
									    b_,
									    yRef.getVector(),
									    izero,
									    sx1Ref.getVector(),
									    izero,
									    gemmR);

									imethod = 2;
									csr_den_kron_mult_method(
									    imethod,
									    transA,
									    transB,
									    a,
									    b_,
									    yRef.getVector(),
									    izero,
									    sx2Ref.getVector(),
									    izero,
									    gemmR);

									imethod = 3;
									csr_den_kron_mult_method(
									    imethod,
									    transA,
									    transB,
									    a,
									    b_,
									    yRef.getVector(),
									    izero,
									    sx3Ref.getVector(),
									    izero,
									    gemmR);

									for (int jx = 0;
									     jx < ncol_X;
									     ++jx) {
										for (int ix = 0;
										     ix < nrow_X;
										     ++ix) {
											auto diff1 = std::abs(
											    x1_(ix,
											        jx)
											    - sx1_(
											        ix,
											        jx));
											auto diff2 = std::abs(
											    x2_(ix,
											        jx)
											    - sx2_(
											        ix,
											        jx));
											auto diff3 = std::abs(
											    x3_(ix,
											        jx)
											    - sx3_(
											        ix,
											        jx));
											auto diffmax = std::max(
											    diff1,
											    std::max(
											        diff2,
											        diff3));
											const double
											    tol
											    = 1.0
											    / (1000.0
											       * 1000.0
											       * 1000.0);
											if (diffmax
											    > tol) {
												INFO(
												    "den_csr: "
												    "itr"
												    "ans"
												    "A "
												    << itransA
												    << "itr"
												       "ans"
												       "B "
												    << itransB
												    << " nr"
												       "ow_"
												       "A "
												    << nrow_A
												    << " nc"
												       "ol_"
												       "A "
												    << ncol_A
												    << " nr"
												       "ow_"
												       "B "
												    << nrow_B
												    << " nco"
												       "l_"
												       "B "
												    << ncol_B
												    << '\n'
												    << "ix "
												    << ix
												    << ", "
												       "jx "
												    << jx
												    << ", "
												       "dif"
												       "f1 "
												    << diff1
												    << ", "
												       "dif"
												       "f2 "
												    << diff2
												    << ", "
												       "dif"
												       "f3 "
												    << diff3
												    << '\n');
												REQUIRE(
												    diffmax
												    <= tol);
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
}

TEST_CASE("kron_mult_test1_double", "[kron][basic]") { run_kron_checks<double>(); }

TEST_CASE("kron_mult_test1_complex", "[kron][basic]") { run_kron_checks<std::complex<double>>(); }

int main(int argc, char* argv[])
{
	Kokkos::initialize(argc, argv);
	int result = Catch::Session().run(argc, argv);
	Kokkos::finalize();
	return result;
}
