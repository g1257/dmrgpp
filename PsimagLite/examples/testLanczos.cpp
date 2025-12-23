#include "CrsMatrix.h"
#include "LanczosSolver.h"
#include "Matrix.h"
#include "PsimagLite.h"

/* This example does Lanczos for a (small) matrix in RAM
 * using PsimagLite's Lanczos solver */
int main(int argc, char* argv[])
{
	PsimagLite::Concurrency concurrency(&argc, &argv, /*nthreads=*/1);
	using ComplexType = std::complex<double>;
	using SolverParametersType = PsimagLite::ParametersForSolver<double>;
	using VectorType = std::vector<ComplexType>;

	/* We fill a dense matrix */
	int n = 8;
	PsimagLite::Matrix<ComplexType> m(n, n);
	// fill m
	for (int i = 0; i < n; ++i)
		m(i, i) = 1.0;
	m(1, 2) = ComplexType(0.0, 0.5);
	m(2, 1) = PsimagLite::conj(m(1, 2));
	m(3, 5) = m(5, 3) = -1.5;

	/* We convert the dense matrix into sparse
	 * Obviously this is unrealistic, we would
	 * start with a CRS matrix from the go.
	 * But here for simplicify of filling I
	 * chose to fill first a dense matrix and then convert,
	 * and also for ease of debugging */
	PsimagLite::CrsMatrix<ComplexType> msparse(m);

	/* These are the parameters that control PsimagLite's Lanczos solver */
	SolverParametersType params;

	/* We store the resulting Lanczos vectors in RAM */
	params.lotaMemory = true;

	/* We create the solver object */
	PsimagLite::
	    LanczosSolver<SolverParametersType, PsimagLite::CrsMatrix<ComplexType>, VectorType>
	        lanczosSolver(msparse, params);

	// This double will contain the lowest eigenvalue
	double e = 0;

	// This will contain the lowest eigenvector
	VectorType z(n, 0.0);

	/* This is the initial vector for Lanczos;
	 * we here set it to random */
	VectorType initial(n);
	PsimagLite::fillRandom(initial);

	// We compute the lowest eigenvalue
	// and eigenvector using Lanczos and starting
	// from vector initial
	// 0 means the lowest
	lanczosSolver.computeOneState(e, z, initial, 0);

	// We print them for comparison and verification
	std::cout << "energy=" << e << "\n";
	std::cout << z << "\n";
}
