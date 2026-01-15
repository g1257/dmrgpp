#ifndef ARNOLDIITERATION_HH
#define ARNOLDIITERATION_HH
#include "Matrix.h"
#include "ParametersForSolver.h"
#include "Sort.h"
#include <algorithm>

namespace PsimagLite {

/* The template parameter MatrixType_ allows the user to pass an object
 * so that the matrix to be decomposed doesn't have to be in memory,
 * only the A*v operation is needed
 */
template <typename MatrixType_> class ArnoldiIteration {
public:

	using MatrixType        = MatrixType_;
	using ComplexOrRealType = typename MatrixType_::value_type;
	using RealType          = typename Real<ComplexOrRealType>::Type;
	using ParamsType        = ParametersForSolver<RealType>;
	using DenseMatrixType   = Matrix<ComplexOrRealType>;
	using VectorType        = std::vector<ComplexOrRealType>;

	ArnoldiIteration(const ParamsType& params)
	    : params_(params)
	{ }

	//---------------------------------------------------------------------------//
	/*!
	 * \brief This function decomposes a into h = q^* a q
	 *
	 * Computes a basis of the (n + 1)-Krylov subspace of the matrix A.
	 * This is the space spanned by the vectors {b, Ab, ..., A^n b},
	 * where b is the initial vector
	 *
	 * \param[out] q The n + 1 vectors of the orthonormal basis
	 *               of the Krylov subspace; each of size m
	 *
	 * \param[out] h  The (n + 1) x n  upper Hessenberg on basis Q.
	 *
	 * \param[in] a The matrix to be decomposed
	 *
	 * \param[in] init_vector The initial vector for Krylov
	 *
	 * \param[in] n The number of Krylov vectors
	 */
	int decompose(std::vector<VectorType>& q,
	              DenseMatrixType&         h,
	              const MatrixType&        a,
	              const VectorType&        init_vector) const
	{
		std::string msg = "ArnoldiIteration::decompose():";
		SizeType    m   = a.rows();
		if (m != a.cols()) {
			throw std::runtime_error(msg + " A not square\n");
		}

		SizeType n = params_.steps;

		h.resize(n, n);
		h.setTo(0.);
		VectorType zero(m);
		q.resize(n, zero);

		double eps       = params_.tolerance;
		int    error_ret = 0;

		// We do not normalize or touch the init_vector
		q[0] = init_vector;

		double     denom = 0; // used as error measure
		VectorType v(m);
		for (SizeType k = 0; k < n; ++k) {
			std::fill(v.begin(), v.end(), 0.);

			a.matrixVectorProduct(v, q[k]); // v += A|q_{k}>

			SizeType total = std::min(n, k + 1);
			for (SizeType j = 0; j < total; ++j) {
				h(j, k) = scalarProduct(q[j], v);
				for (SizeType jj = 0; jj < m; ++jj) {
					v[jj] -= h(j, k) * q[j][jj];
				}
			}

			// Make h square so we can diagonalize
			if (k + 1 == n) {
				break;
			}

			h(k + 1, k) = PsimagLite::norm(v);
			denom       = PsimagLite::real(h(k + 1, k));
			if (eps > 0) { // User wants early exit if possible (recommended)
				if (denom > eps) {
					double              factor = 1. / denom;
					q[k + 1] <= factor* v;
				} else {
					q.resize(k + 1);
					h.resize(k + 1, k + 1);
					break;
				}
			} else { // Do all iterations, no early exit
				bool is_zero = (denom == 0.);
				if (is_zero && error_ret == 0) {
					error_ret = k; // k is assured to be greater than zero
				}

				double              factor = (is_zero) ? 1. : 1. / denom;
				q[k + 1] <= factor* v;
			}
		}

		std::cout << "ArnolidIteration: done after " << q.size() << " steps;";
		std::cout << " error=" << denom << "\n";
		return error_ret;
	}

	//---------------------------------------------------------------------------//
	/*!
	 * \brief Diagonalize the Hessenberg H
	 *
	 * Note that the Hessenberg H is not necessarily Hermitian so
	 * that the eigs_of_h are always considered vector<complex> even
	 * if h is real.
	 *
	 * \param[out] eigs_of_h             The eigenvalues of H
	 * \param[out] h_right_eigenvectors  The right eigenvectors of H
	 * \param[in]  h                     The hessenberg matrix
	 */
	void diagHessenberg(std::vector<std::complex<RealType>>&         eigs_of_h,
	                    PsimagLite::Matrix<ComplexOrRealType>&       h_right_eigenvectors,
	                    const PsimagLite::Matrix<ComplexOrRealType>& h) const
	{
		SizeType n_arnoldi = h.rows();
		if (n_arnoldi != h.cols()) {
			throw RuntimeError("diagHessenberg(): Hessenberg not square\n");
		}

		eigs_of_h.resize(n_arnoldi);
		h_right_eigenvectors.resize(n_arnoldi, n_arnoldi);
		PsimagLite::Matrix<ComplexOrRealType> vl_unused(n_arnoldi, n_arnoldi);

		PsimagLite::Matrix<ComplexOrRealType> h_copy(h);
		geev('N', 'V', h_copy, eigs_of_h, vl_unused, h_right_eigenvectors);
	}

	//---------------------------------------------------------------------------//
	/*!
	 * \brief Get Arnoldi eigenvector number col from the Arnoldi decomposition
	 *
	 * \param[out] ev  The eigenvector
	 * \param[in]  s   The matrix that diagonalizes the Hessenberg
	 * \param[in]  q   The Q matrix
	 * \param[in]  col The index of the eigenvector wanted
	 */
	template <typename ComplexOrRealType>
	void getEigenvector(std::vector<ComplexOrRealType>&                    ev,
	                    const PsimagLite::Matrix<ComplexOrRealType>&       s,
	                    const std::vector<std::vector<ComplexOrRealType>>& q,
	                    SizeType                                           col) const
	{
		if (q.empty()) {
			throw RuntimeError("getArnolidEigenvector(): q.size==0\n");
		}

		unsigned int nbig   = q[0].size();
		unsigned int nsmall = q.size();

		if (nsmall != s.rows() || nsmall != s.cols()) {
			throw RuntimeError("getArnolidEigenvector(): q.size != s.rows or s.cols\n");
		}

		if (col >= s.cols()) {
			throw RuntimeError("getArnolidEigenvector(): col >= s.cols\n");
		}

		for (SizeType i = 0; i < nbig; ++i) {
			ComplexOrRealType sum = 0;
			for (SizeType j = 0; j < nsmall; ++j) {
				sum += q[j][i] * s(j, col);
			}

			ev[i] = sum;
		}
	}

private:

	ParamsType params_;
};
}
#endif // ARNOLDIITERATION_HH
