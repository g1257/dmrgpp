#include "ArnoldiIteration.hh"
#include "CrsMatrix.h"
#include "Matrix.h"
#include "PsimagLite.h"
#include "Random48.h"
#include <cassert>

using ComplexType          = std::complex<double>;
using ComplexOrRealType    = ComplexType;
using SolverParametersType = PsimagLite::ParametersForSolver<double>;
using VectorType           = std::vector<ComplexOrRealType>;
using RandomType           = PsimagLite::Random48<double>;

/* This example does Arnoldi iteration
 * using PsimagLite's ArnoldiSaI solver */
void test_arnoldi()
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

	// Diag the matrix for testing (FIXME: move to a function)
	PsimagLite::Matrix<ComplexOrRealType> m_copy(m); // m_copy will be overwritten by geev
	std::vector<ComplexType>              eigenvalues(n);
	PsimagLite::Matrix<ComplexOrRealType> vl_unused(n, n);
	PsimagLite::Matrix<ComplexOrRealType> right_eigenvectors(n, n);
	geev('N', 'V', m_copy, eigenvalues, vl_unused, right_eigenvectors);

	// Matrix eigenvalues are these
	std::cout << "Actual eigs\n";
	std::cout << eigenvalues << "\n\n";

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
	std::cerr << "Arnolid return code = " << error << "\n";

	// Now we diag h
	std::vector<ComplexOrRealType>        eigs_of_h;
	PsimagLite::Matrix<ComplexOrRealType> h_right_eigenvectors;
	arnoldi.diagHessenberg(eigs_of_h, h_right_eigenvectors, h);
	std::cout << "Eigs of h\n";
	std::cout << eigs_of_h << "\n\n";

	std::vector<ComplexOrRealType> l_eigenvector_estimate(n);
	arnoldi.getEigenvector(l_eigenvector_estimate, h_right_eigenvectors, q, 0);
	// makeItReal(l_eigenvector_estimate);
	std::cout << "\nEstimated eigenvector \n";
	std::cout << l_eigenvector_estimate << "\n";
}

int main(int argc, char* argv[])
{
	PsimagLite::Concurrency concurrency(&argc, &argv, /*nthreads=*/1);
	test_arnoldi();
}
