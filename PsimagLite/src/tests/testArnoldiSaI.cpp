#include "ArnoldiSaI.hh"
#include "CrsMatrix.h"
#include "Matrix.h"
#include "PsimagLite.h"
#include <cassert>
#include <catch2/catch_approx.hpp>
#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>
#include <random>

using ComplexType          = std::complex<double>;
using ComplexOrRealType    = ComplexType;
using SolverParametersType = PsimagLite::ParametersForSolver<double>;
using VectorType           = std::vector<ComplexOrRealType>;
using RandomType           = PsimagLite::Random48<double>;

//---------------------------------------------------------------------------//
/*!
 * \brief Finds the index of the element with lowest complex norm
 *
 * \param[in] eigs A vector of complex numbers
 *
 * \returns the index of the element with lowest complex norm
 */
unsigned int findTheOne(const std::vector<std::complex<double>>& eigs)
{
	auto min_element_it
	    = std::min_element(eigs.begin(),
	                       eigs.end(),
	                       [](const std::complex<double>& a, const std::complex<double>& b)
	                       { return std::norm(a) < std::norm(b); });

	assert(min_element_it != eigs.end());
	return min_element_it - eigs.begin();
}

/* This example does Arnoldi shift-and-invert for a (small) matrix in RAM
 * using PsimagLite's ArnoldiSaI solver */
TEST_CASE("Full Arnoldi shift-and-invert of a random matrix", "[ArnoldiSaI]")
{
	/* We fill a dense matrix */
	int                                   n   = 64;
	double                                max = 10.;
	double                                min = -4.;
	PsimagLite::Matrix<ComplexOrRealType> m(n, n);
	// fill m
	// A fixed seed value for reproducibility
	constexpr unsigned int SEED = 12345;

	// 1. Seed the random number engine with the fixed value
	std::mt19937                           rng(SEED);
	std::uniform_real_distribution<double> dist(min, max);

	for (int i = 0; i < n; ++i) {
		for (int j = 0; j < n; ++j) {
			m(i, j) = dist(rng);
		}

		m(i, i) = 0.;
	}

	// Diag the matrix for testing
	PsimagLite::Matrix<ComplexOrRealType> m_copy(m); // m_copy will be overwritten by geev
	std::vector<ComplexType>              eigenvalues(n);
	PsimagLite::Matrix<ComplexOrRealType> vl_unused(n, n);
	PsimagLite::Matrix<ComplexOrRealType> right_eigenvectors(n, n);
	geev('N', 'V', m_copy, eigenvalues, vl_unused, right_eigenvectors);

	REQUIRE(n == eigenvalues.size());
	REQUIRE_THAT(0., Catch::Matchers::WithinAbs(std::abs(std::imag(eigenvalues[0])), 1e-5));

	unsigned int the_one = findTheOne(eigenvalues);
	assert(the_one < eigenvalues.size());

	// Now we do Arnoldi Shift-and-Invert
	double sigma = 0.1; // small value here

	// shift first
	m_copy = m;
	for (int i = 0; i < n; ++i) {
		m_copy(i, i) -= sigma;
	}

	// invert now
	inverse(m_copy);

	// Eigenvalues of the inverse
	PsimagLite::Matrix<ComplexOrRealType> matrix_inverse(m_copy);
	std::vector<ComplexType>              inverse_eigs(n);
	geev('N', 'V', m_copy, inverse_eigs, vl_unused, right_eigenvectors);

	// Predict for comparison the eigenvalues of the sai
	ComplexType predicted_eig = 1.0 / (eigenvalues[the_one] - sigma);
	int         index_inverse = n - 1 - the_one;
	REQUIRE(index_inverse < inverse_eigs.size());
	ComplexType value_inverse = inverse_eigs[index_inverse];
	CHECK(std::real(predicted_eig) == Catch::Approx(std::real(value_inverse)));
	CHECK_THAT(std::abs(std::imag(predicted_eig)),
	           Catch::Matchers::WithinAbs(std::abs(std::imag(value_inverse)), 1e-5));

	// Finally do ArnoldiSaI
	/* We convert the dense matrix into sparse
	 * Obviously this is unrealistic, we would
	 * start with a CRS matrix from the go.
	 * But here for simplicify of filling I
	 * chose to fill first a dense matrix and then convert,
	 * and also for ease of debugging */
	PsimagLite::CrsMatrix<ComplexOrRealType> msparse(m);

	/* These are the parameters that control PsimagLite's Lanczos solver */
	SolverParametersType params;

	/* We store the resulting Lanczos vectors in RAM */
	params.lotaMemory = true;
	params.tolerance  = 1e-8;
	params.steps      = 50;

	/* We create the solver object */
	PsimagLite::ArnoldiSaI<PsimagLite::CrsMatrix<ComplexOrRealType>> arnoldi_sai(
	    msparse, params, sigma);

	/* This is the initial vector for Lanczos;
	 * we here set it to random */
	VectorType initial(n);
	PsimagLite::fillRandom(initial);

	std::vector<ComplexOrRealType> eigenvector;
	double                         eigenvalue  = 0;
	constexpr SizeType             bogus_index = 0;
	arnoldi_sai.computeOneState(eigenvalue, eigenvector, initial, bogus_index);
	CHECK(std::real(eigenvalue) == Catch::Approx(std::real(eigenvalues[the_one])));
}
