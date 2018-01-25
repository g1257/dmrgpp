#ifndef MATRIXDENSEORSPARSE_H
#define MATRIXDENSEORSPARSE_H
#include "Vector.h"
#include "KronUtilWrapper.h"
#include "Matrix.h"

namespace Dmrg {

template<typename SparseMatrixType>
class MatrixDenseOrSparse {

public:

	typedef typename SparseMatrixType::value_type ComplexOrRealType;
	typedef typename PsimagLite::Real<ComplexOrRealType>::Type RealType;
	typedef typename PsimagLite::Vector<ComplexOrRealType>::Type VectorType;
	typedef PsimagLite::Vector<int>::Type VectorIntType;

	explicit MatrixDenseOrSparse(const SparseMatrixType& sparse,
	                             const RealType& threshold)
	    : isDense_(sparse.nonZero() > static_cast<int>(threshold*sparse.rows()*sparse.cols())),
	      sparseMatrix_(sparse)
	{
		sparseMatrix_.checkValidity();

		if (isDense_) // A(i,j) at  val[ (i) + (j)*nrow ]
			crsMatrixToFullMatrix(denseMatrix_, sparse);
	}

	bool isDense() const { return isDense_; }

	SizeType rows() const
	{
		return sparseMatrix_.rows();
	}

	SizeType cols() const
	{
		return sparseMatrix_.cols();
	}

	const PsimagLite::Matrix<ComplexOrRealType>& dense() const
	{
		if (!isDense_)
			throw PsimagLite::RuntimeError("FATAL: Matrix isn't dense\n");
		return denseMatrix_;
	}

	const SparseMatrixType& sparse() const
	{
		sparseMatrix_.checkValidity();
		return sparseMatrix_;
	}

	bool isZero() const
	{
		return (isDense_) ? false : (sparseMatrix_.nonZero() == 0);
	}

	SparseMatrixType toSparse() const
	{
		return (isDense_) ? SparseMatrixType(denseMatrix_) : sparse();
	}

private:

	bool isDense_;
	const PsimagLite::CrsMatrix<ComplexOrRealType> sparseMatrix_;
	PsimagLite::Matrix<ComplexOrRealType> denseMatrix_;
}; // class MatrixDenseOrSparse

template<typename SparseMatrixType>
void kronMult(typename PsimagLite::Vector<typename SparseMatrixType::value_type>::Type& xout,
              SizeType offsetX,
              const typename PsimagLite::Vector<typename SparseMatrixType::value_type>::Type& yin,
              SizeType offsetY,
              char transA,
              char transB,
              const MatrixDenseOrSparse<SparseMatrixType>& A,
              const MatrixDenseOrSparse<SparseMatrixType>& B)
{
	const bool isDenseA = A.isDense();
	const bool isDenseB = B.isDense();

	if (isDenseA) {
		if (isDenseB) {
			den_kron_mult(transA,
			              transB,
			              A.dense(),
			              B.dense(),
			              yin,
			              offsetY,
			              xout,
			              offsetX);
		} else  {
			// B is sparse
			den_csr_kron_mult(transA,
			                  transB,
			                  A.dense(),
			                  B.sparse(),
			                  yin,
				              offsetY,
				              xout,
				              offsetX);
		}
	} else {
		// A is sparse
		if (isDenseB) {
			csr_den_kron_mult(transA,
			                  transB,
			                  A.sparse(),
			                  B.dense(),
			                  yin,
				              offsetY,
				              xout,
				              offsetX);
		} else {
			// B is sparse
			csr_kron_mult(transA,
			              transB,
			              A.sparse(),
			              B.sparse(),
			              yin,
			              offsetY,
			              xout,
			              offsetX);
		};
	};
} // kron_mult

} // namespace Dmrg
#endif // MATRIXDENSEORSPARSE_H
