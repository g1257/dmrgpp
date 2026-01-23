#include "ArnoldiIteration.hh"
#include "CrsMatrix.h"
#include "Matrix.h"
#include "PsimagLite.h"
#include <cassert>
#include <catch2/catch_approx.hpp>
#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>
#include <random>

namespace PsimagLite {

using ComplexType          = std::complex<double>;
using ComplexOrRealType    = ComplexType;
using SolverParametersType = PsimagLite::ParametersForSolver<double>;
using VectorType           = std::vector<ComplexOrRealType>;
using RandomType           = PsimagLite::Random48<double>;

/* This example does Arnoldi iteration
 * using PsimagLite's ArnoldiSaI solver */
TEST_CASE("Full ArnoldiIteration of a random matrix", "[ArnoldiIteration]")
{
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
	}

	// Diag the matrix for testing
	PsimagLite::Matrix<ComplexOrRealType> m_copy(m); // m_copy will be overwritten by geev
	std::vector<ComplexType>              eigenvalues(n);
	PsimagLite::Matrix<ComplexOrRealType> vl_unused(n, n);
	PsimagLite::Matrix<ComplexOrRealType> right_eigenvectors(n, n);

	geev('N', 'V', m_copy, eigenvalues, vl_unused, right_eigenvectors);
	REQUIRE(eigenvalues.size() == n);
	ComplexType ref_eig = eigenvalues[0];
	CHECK_THAT(std::abs(std::imag(ref_eig)), Catch::Matchers::WithinAbs(0, 1e-5));

	// Finally do Arnoldi
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
	PsimagLite::ArnoldiIteration<PsimagLite::CrsMatrix<ComplexOrRealType>> arnoldi(params);

	/* This is the initial vector for Lanczos;
	 * we here set it to random */
	VectorType initial(n);
	PsimagLite::fillRandom(initial);

	std::vector<VectorType>               q;
	PsimagLite::Matrix<ComplexOrRealType> h;
	int                                   error = arnoldi.decompose(q, h, msparse, initial);
	CHECK(0 == error);

	// Now we diag h
	std::vector<ComplexOrRealType>        eigs_of_h;
	PsimagLite::Matrix<ComplexOrRealType> h_right_eigenvectors;
	arnoldi.diagHessenberg(eigs_of_h, h_right_eigenvectors, h);

	REQUIRE(eigs_of_h.size() > 0);
	CHECK(std::real(ref_eig) == Catch::Approx(std::real(eigs_of_h[0])));
	CHECK_THAT(0., Catch::Matchers::WithinAbs(std::abs(std::imag(ref_eig)), 1e-5));
	std::vector<ComplexOrRealType> l_eigenvector_estimate(n);
	arnoldi.getEigenvector(l_eigenvector_estimate, h_right_eigenvectors, q, 0);
}
}
