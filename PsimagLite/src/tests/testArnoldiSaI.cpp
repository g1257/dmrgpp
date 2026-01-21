#include "ArnoldiSaI.hh"
#include "CrsMatrix.h"
#include "Matrix.h"
#include "PsimagLite.h"
#include "Random48.h"
#include <cassert>
#include <catch2/catch_approx.hpp>
#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

using ComplexType          = std::complex<double>;
using ComplexOrRealType    = ComplexType;
using SolverParametersType = PsimagLite::ParametersForSolver<double>;
using VectorType           = std::vector<ComplexOrRealType>;
using RandomType           = PsimagLite::Random48<double>;

//---------------------------------------------------------------------------//
/*!
 * \brief Lowest element that is also real and positive
 *
 * \param[in] eigs A vector of complex numbers
 *
 * \returns the lowest element that is also real and positive
 */
unsigned int findTheOne(const std::vector<std::complex<double>>& eigs)
{
	constexpr double          eps = 1e-8; // to determine real
	unsigned int              n   = eigs.size();
	std::vector<double>       candidates(n); // preallocate
	std::vector<unsigned int> indices(n); // preallocate
	unsigned int              j = 0;
	for (unsigned int i = 0; i < n; ++i) {
		const std::complex<double>& val = eigs[i];
		// real?
		if (std::abs(std::imag(val)) > eps) {
			continue;
		}

		double x = std::real(val);
		// postive?
		if (x < 0) {
			continue;
		}

		candidates[j] = x;
		indices[j++]  = i;
	}

	candidates.resize(j);
	indices.resize(j);
	PsimagLite::Sort<std::vector<double>> sort;
	std::vector<long unsigned int>        iperm(j);
	sort.sort(candidates, iperm);

	return indices[iperm[0]];
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
	constexpr long int seed = 1234;
	RandomType         rng(seed);
	for (int i = 0; i < n; ++i) {
		for (int j = 0; j < n; ++j) {
			m(i, j) = rng.random() * max + min;
		}
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
	ComplexType predicted_eig = 1.0 / (eigenvalues[0] - sigma);
	int         index_inverse = n - 1;
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
