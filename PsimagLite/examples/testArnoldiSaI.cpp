#include "ArnoldiSaI.hh"
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

//---------------------------------------------------------------------------//
/*!
 * \brief Print matrix column
 *
 * \param[out] os     The stream to print to
 * \param[in]  mv     The matrix of complex or real values
 * \param[in]  index  The column index to print
 */
template <typename ComplexOrRealType>
void printEigenvector(std::ostream&                                os,
                      const PsimagLite::Matrix<ComplexOrRealType>& mv,
                      unsigned int                                 index)
{
	os << "Eigenvector for eigenvalue index = " << index << "\n";
	SizeType n = mv.rows();
	for (SizeType i = 0; i < n; ++i) {
		os << mv(i, index) << " ";
	}

	os << "\n";
}

/* This example does Arnoldi shift-and-invert for a (small) matrix in RAM
 * using PsimagLite's ArnoldiSaI solver */
void test_arnoldi_sai()
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

		// m(i, i) = 0;
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

	unsigned int the_one = findTheOne(eigenvalues);
	assert(the_one < eigenvalues.size());
	std::cout << "Lowest eig of L = " << the_one << " " << eigenvalues[the_one] << "\n";
	printEigenvector(std::cout, right_eigenvectors, the_one);

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
	std::cout << "Inverse eigs\n";
	std::cout << inverse_eigs << "\n\n";

	// Predict for comparison the eigenvalues of the sai
	std::vector<ComplexType> predicted_eigs(n);
	for (int i = 0; i < n; ++i) {
		predicted_eigs[i] = 1.0 / (eigenvalues[i] - sigma);
	}

	std::cout << "Predicted eigs\n";
	std::cout << predicted_eigs << "\n\n";

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
	std::cout << "Eigenvalue = " << eigenvalue << "\n";

	std::cout << "\nEstimated eigenvector \n";
	std::cout << eigenvector << "\n";
}

int main(int argc, char* argv[])
{
	PsimagLite::Concurrency concurrency(&argc, &argv, /*nthreads=*/1);
	test_arnoldi_sai();
}
