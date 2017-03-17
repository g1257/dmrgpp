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
	typedef typename PsimagLite::Vector<ComplexOrRealType>::Type VectorType;
	typedef PsimagLite::Vector<int>::Type VectorIntType;

	// 50% cutoff == 0.5 here for sparse/dense
	explicit MatrixDenseOrSparse(const SparseMatrixType& sparse)
	    : isDense_(true), //sparse.nonZero() > static_cast<int>(0.5*sparse.row()*sparse.col())),
	      sparseMatrix_(0),
	      denseMatrix_(sparse.row(), sparse.col())
	{
		if (!isDense_) {
			sparseMatrix_ = &sparse;
		} else {
			// A(i,j) at  val[ (i) + (j)*nrow ]
			crsMatrixToFullMatrix(denseMatrix_, sparse);
		}
	}

	bool isDense() const { return isDense_; }

	SizeType rows() const { return denseMatrix_.n_row(); }

	SizeType cols() const { return denseMatrix_.n_col(); }

	const PsimagLite::Matrix<ComplexOrRealType>& dense() const
	{
		return denseMatrix_;
	}

	const SparseMatrixType& sparse() const
	{
		assert(!sparseMatrix_);
		return *sparseMatrix_;
	}

	bool isZero() const
	{
		return (isDense_) ? false : (sparseMatrix_.nonZero() == 0);
	}

	SparseMatrixType toSparse() const
	{
		if (isDense_) {
			return SparseMatrixType(denseMatrix_);
		}

		return sparse();
	}

private:

	bool isDense_;
	const PsimagLite::CrsMatrix<ComplexOrRealType>* sparseMatrix_;
	PsimagLite::Matrix<ComplexOrRealType> denseMatrix_;
}; // class MatrixDenseOrSparse

template<typename SparseMatrixType>
void kronMult(typename SparseMatrixType::value_type* xout,
              const typename SparseMatrixType::value_type* yin,
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
			              xout);
		} else  {
			// B is sparse
			den_csr_kron_mult(transA,
			                  transB,
			                  A.dense(),
			                  B.sparse(),
			                  yin,
			                  xout);
		}
	} else {
		// A is sparse
		if (isDenseB) {
			csr_den_kron_mult(transA,
			                  transB,
			                  A.sparse(),
			                  B.dense(),
			                  yin,
			                  xout);
		} else {
			// B is sparse
			csr_kron_mult(transA,
			              transB,
			              A.sparse(),
			              B.sparse(),
			              yin,
			              xout);
		};
	};
} // kron_mult

} // namespace Dmrg
#endif // MATRIXDENSEORSPARSE_H
