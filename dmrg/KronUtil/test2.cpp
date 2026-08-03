#include "KronUtil.h"
#include "util.h"

#include <Kokkos_Core.hpp>
#define CATCH_CONFIG_RUNNER
#include <catch2/catch_session.hpp>
#include <catch2/catch_test_macros.hpp>
#include <complex>

template <typename ComplexOrRealType> void run_kron_submatrix_checks()
{
	using RealType = typename PsimagLite::Real<ComplexOrRealType>::Type;

	for (int thresholdB_idx = 0; thresholdB_idx <= 11; ++thresholdB_idx) {
		double thresholdB = .1 * thresholdB_idx;
		for (int thresholdA_idx = 0; thresholdA_idx <= 11; ++thresholdA_idx) {
			double thresholdA = .1 * thresholdA_idx;
			for (int ncol_A = 1; ncol_A <= 10; ncol_A += 3) {
				for (int nrow_A = 1; nrow_A <= 10; nrow_A += 3) {
					for (int ncol_B = 1; ncol_B <= 10; ncol_B += 3) {
						for (int nrow_B = 1; nrow_B <= 10; nrow_B += 3) {
							PsimagLite::Matrix<ComplexOrRealType> a_(
							    nrow_A, ncol_A);
							PsimagLite::Matrix<ComplexOrRealType> b_(
							    nrow_B, ncol_B);

							if ((thresholdA == 0)
							    && (nrow_A == ncol_A)) {
								/*
								 * ------------------------------------
								 * special case to test identity
								 * matrix
								 * ------------------------------------
								 */
								den_eye(nrow_A, ncol_A, a_);
								REQUIRE(den_is_eye(a_));
								PsimagLite::CrsMatrix<
								    ComplexOrRealType>
								    a(a_);
								REQUIRE(csr_is_eye(a));
							} else {
								den_gen_matrix(
								    nrow_A, ncol_A, thresholdA, a_);

								PsimagLite::CrsMatrix<
								    ComplexOrRealType>
								    a(a_);
								REQUIRE(den_is_eye(a_)
								        == csr_is_eye(a));
								REQUIRE(den_is_zeros(a_)
								        == csr_is_zeros(a));
							}

							if ((thresholdB == 0)
							    && (nrow_B == ncol_B)) {
								den_eye(nrow_B, ncol_B, b_);
								REQUIRE(den_is_eye(b_));
								PsimagLite::CrsMatrix<
								    ComplexOrRealType>
								    b(b_);
								REQUIRE(csr_is_eye(b));
							} else {
								den_gen_matrix(
								    nrow_B, ncol_B, thresholdB, b_);

								PsimagLite::CrsMatrix<
								    ComplexOrRealType>
								    b(b_);
								REQUIRE(den_is_eye(b_)
								        == csr_is_eye(b));
								REQUIRE(den_is_zeros(b_)
								        == csr_is_zeros(b));
							}
							/*
							 * -----------------------------------
							 * generate compressed row format from
							 * dense matrices A, B
							 * -----------------------------------
							 */
							PsimagLite::CrsMatrix<ComplexOrRealType> a(
							    a_);
							REQUIRE(den_is_eye(a_) == csr_is_eye(a));

							PsimagLite::CrsMatrix<ComplexOrRealType> b(
							    b_);
							REQUIRE(den_is_eye(a_) == csr_is_eye(a));

							/*
							 * -----------------------------
							 * explicitly form C = kron(A,B)
							 * -----------------------------
							 */

							const int nrow_C = nrow_A * nrow_B;
							const int ncol_C = ncol_A * ncol_B;
							if (ncol_C < 2)
								continue;
							PsimagLite::Matrix<ComplexOrRealType> c_(
							    nrow_C, ncol_C);

							den_kron_form(nrow_A,
							              ncol_A,
							              a_,
							              nrow_B,
							              ncol_B,
							              b_,
							              c_);

							/*
							 * ---------------------------------------
							 * generate compressed sparse version of C
							 * ---------------------------------------
							 */
							PsimagLite::CrsMatrix<ComplexOrRealType> c(
							    c_);

							PsimagLite::Vector<int>::Type rindex(
							    nrow_C);
							PsimagLite::Vector<int>::Type cindex(
							    ncol_C);

							/*
							 * --------------------
							 * extract  even rows
							 * extract odd columns
							 * --------------------
							 */
							int nrindex = 0;
							for (int ic = 0; ic < nrow_C; ic += 2) {
								rindex[nrindex++] = ic;
							}
							int ncindex = 0;
							for (int jc = 1; jc < ncol_C; jc += 2) {
								cindex[ncindex++] = jc;
							}

							/*
							 * -------------------------------
							 * extract submatrix from C into D
							 * -------------------------------
							 */

							int nrow_D = nrindex;
							int ncol_D = ncindex;
							PsimagLite::Matrix<ComplexOrRealType> d_(
							    nrow_D, ncol_D);

							den_submatrix(nrow_C,
							              ncol_C,
							              c_,
							              nrindex,
							              ncindex,
							              rindex,
							              cindex,
							              d_);
							/*
							 * -------------------------------------
							 * generate submatrix from sparse version C
							 * -------------------------------------
							 */

							const int max_nnz_D = 1 + den_nnz(d_);

							PsimagLite::CrsMatrix<ComplexOrRealType> sd(
							    d_.n_row(), d_.n_col());

							csr_submatrix(c,

							              nrindex,
							              ncindex,
							              max_nnz_D,
							              rindex,
							              cindex,

							              sd);

							/*
							 * -----------------------
							 * convert to dense matrix
							 * -----------------------
							 */
							PsimagLite::Matrix<ComplexOrRealType> dd_(
							    nrow_D, ncol_D);
							crsMatrixToFullMatrix(dd_, sd);

							/*
							 * --------------------------------------
							 * check both matrices should be the same
							 * --------------------------------------
							 */

							for (int jd = 0; jd < ncol_D; ++jd) {
								for (int id = 0; id < nrow_D;
								     ++id) {
									RealType diff = std::abs(
									    dd_(id, jd)
									    - d_(id, jd));
									const RealType tol = 1.0
									    / (1000.0 * 1000.0);
									if (diff > tol) {
										INFO("nrow_D "
										     << nrow_D
										     << " ncol_D "
										     << ncol_D
										     << " DD(" << id
										     << "," << jd
										     << ")"
										     << dd_(id, jd)
										     << " D(" << id
										     << "," << jd
										     << ")"
										     << d_(id, jd)
										     << " di"
										        "ff "
										     << diff);
										REQUIRE(diff
										        <= tol);
									}
								}
							}

							/*
							 * --------------------------------
							 * form E = C(rindex(:),cindex(:))
							 * without forming C
							 * E should be the same as matrix D
							 * --------------------------------
							 */
							int nrow_E = nrindex;
							int ncol_E = ncindex;
							PsimagLite::Matrix<ComplexOrRealType> e_(
							    nrow_E, ncol_E);

							den_kron_submatrix(nrow_A,
							                   ncol_A,
							                   a_,
							                   nrow_B,
							                   ncol_B,
							                   b_,
							                   nrindex,
							                   ncindex,
							                   rindex,
							                   cindex,
							                   e_);

							/*
							 * --------------------------
							 * check E and D are the same
							 * --------------------------
							 */
							for (int je = 0; je < ncol_E; ++je) {
								for (int ie = 0; ie < nrow_E;
								     ++ie) {
									ComplexOrRealType eij
									    = e_(ie, je);
									ComplexOrRealType dij
									    = d_(ie, je);

									if (eij != dij) {
										INFO("nrow_A="
										     << nrow_A
										     << " ncol_A="
										     << ncol_A
										     << " nrow_B="
										     << nrow_B
										     << " ncol_B="
										     << ncol_B
										     << "  nrindex "
										     << nrindex
										     << " ncindex "
										     << ncindex
										     << '\n'
										     << " ie " << ie
										     << " je " << je
										     << " eij "
										     << eij
										     << " dij "
										     << dij);
										REQUIRE(eij == dij);
									}
								}
							}

							/*
							 * ------------------------------
							 * form sparse version of E from
							 * sparse version of A, B
							 * ------------------------------
							 */

							PsimagLite::CrsMatrix<ComplexOrRealType> e(
							    nrow_E, ncol_E);
							const int max_nnz_E = 1 + den_nnz(e_);
							csr_kron_submatrix(a,
							                   b,
							                   nrindex,
							                   ncindex,
							                   max_nnz_E,
							                   rindex,
							                   cindex,
							                   e);

							/*
							 * ---------------------------------
							 * convert from sparse back to dense
							 * ---------------------------------
							 */
							PsimagLite::Matrix<ComplexOrRealType> se_(
							    nrow_E, ncol_E);

							crsMatrixToFullMatrix(se_, e);

							/*
							 * -----------------------------------------
							 * check E(ie,je) and SE(ie,je) are the same
							 * -----------------------------------------
							 */

							for (int je = 0; je < ncol_E; ++je) {
								for (int ie = 0; ie < nrow_E;
								     ++ie) {
									RealType diff = std::abs(
									    e_(ie, je)
									    - se_(ie, je));
									const RealType tol = 1.0
									    / (1000.0 * 1000.0);
									if (diff > tol) {
										INFO(
										    "nrow_E="
										    << nrow_E
										    << " ncol_E="
										    << ncol_E

										    << " E(" << ie
										    << "," << je
										    << " "
										    << e_(ie, je)
										    << " SE(" << ie
										    << "," << je
										    << ")"
										    << se_(ie, je));
										REQUIRE(diff
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

TEST_CASE("kron_submatrix_test2_double", "[kron][submatrix]")
{
	run_kron_submatrix_checks<double>();
}
TEST_CASE("kron_submatrix_test2_complex", "[kron][submatrix]")
{
	run_kron_submatrix_checks<std::complex<double>>();
}

int main(int argc, char* argv[])
{
	Kokkos::initialize(argc, argv);
	int result = Catch::Session().run(argc, argv);
	Kokkos::finalize();
	return result;
}
