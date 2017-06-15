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

	explicit MatrixDenseOrSparse(const VectorType& v,
	                             SizeType rows,
	                             SizeType cols,
	                             const RealType& threshold)
	    : rows_(rows), cols_(cols), denseMatrix_(v, rows, cols)
	{		
		isDense_ = (denseMatrix_.nonZeros() > threshold*rows*cols);
		if (isDense_) return;
		fullMatrixToCrsMatrix(sparseMatrix_, denseMatrix_);
	}

	bool isDense() const { return isDense_; }

	SizeType rows() const
	{
		return rows_;
	}

	SizeType cols() const
	{
		return cols_;
	}

	const PsimagLite::Matrix<ComplexOrRealType>& dense() const
	{
		return denseMatrix_;
	}

	const SparseMatrixType& sparse() const
	{
		if (isDense_)
			err("MatrixDenseOrSparse::sparse() cannot be called when isDense\n");

		return sparseMatrix_;
	}

	bool isZero() const
	{
		return (isDense_) ? false : (sparseMatrix_.nonZero() == 0);
	}

private:

	SizeType rows_;
	SizeType cols_;
	bool isDense_;
	PsimagLite::CrsMatrix<ComplexOrRealType> sparseMatrix_;
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
