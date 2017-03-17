#ifndef MATRIXDENSEORSPARSE_H
#define MATRIXDENSEORSPARSE_H
#include "Vector.h"
#include "KronUtilWrapper.h"

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
	      denseMatrix_(sparse.row(), sparse.col())
	{
		SizeType rows = sparse.row();
		if (!isDense_) {
			rowptr_.resize(rows + 1, 0);
			for (SizeType k = 0; k < rows + 1; ++k)
				rowptr_[k] = sparse.getRowPtr(k);

			SizeType nz = sparse.nonZero();
			if (nz == 0) return;
			colind_.resize(nz, 0);
			values_.resize(nz, 0.0);
			for (SizeType k = 0; k < nz; ++k) {
				colind_[k] = sparse.getCol(k);
				values_[k] = sparse.getValue(k);
			}

		} else {
			// A(i,j) at  val[ (i) + (j)*nrow ]
			for (SizeType i = 0; i < rows; ++i) {
				for (int k = sparse.getRowPtr(i); k < sparse.getRowPtr(i+1); ++k)
					denseMatrix_(i, sparse.getCol(k)) = sparse.getValue(k);
			}
		}
	}

	bool isDense() const { return isDense_; }

	SizeType rows() const { return denseMatrix_.n_row(); }

	SizeType cols() const { return denseMatrix_.n_col(); }

	const PsimagLite::Matrix<ComplexOrRealType>& toDense() const
	{
		return denseMatrix_;
	}

	const VectorType& values() const
	{
		assert(!isDense_);
		return values_;
	}

	const VectorIntType& colind() const
	{
		assert(!isDense());
		return colind_;
	}

	const VectorIntType& rowptr() const
	{
		assert(!isDense());
		return rowptr_;
	}

	bool isZero() const
	{
		return (isDense_) ? false : (values_.size() == 0);
	}

	SparseMatrixType toSparse() const
	{
		if (isDense_) {
			return SparseMatrixType(denseMatrix_);
		}

		SparseMatrixType m;
		SizeType rows = this->rows();
		m.resize(rows, cols(), colind_.size());
		SizeType counter = 0;
		for (SizeType i = 0; i < rows; ++i) {
			m.setRow(i, counter);
			for (int k = rowptr_[i]; k < rowptr_[i+1]; ++k)
				m.setCol(counter++,colind_[k]);
		}

		m.setRow(rows, counter);
		m.checkValidity();
		return m;
	}

private:

	bool isDense_;
	VectorIntType rowptr_;
	VectorIntType colind_;
	VectorType values_;
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
	SizeType nrowA = A.rows();
	SizeType ncolA = A.cols();
	SizeType nrowB = B.rows();
	SizeType ncolB = B.cols();
	const bool isDenseA = A.isDense();
	const bool isDenseB = B.isDense();

	if (isDenseA) {
		if (isDenseB) {
			den_kron_mult(transA,
			              transB,
			              A.toDense(),
			              B.toDense(),
			              yin,
			              xout);
		} else  {
			// B is sparse
			den_csr_kron_mult(transA,
			                  transB,
			                  A.toDense(),
			                  nrowB,
			                  ncolB,
			                  B.rowptr(),
			                  B.colind(),
			                  B.values(),
			                  yin,
			                  xout);
		}
	} else {
		// A is sparse
		if (isDenseB) {
			csr_den_kron_mult(transA,
			                  transB,
			                  nrowA,
			                  ncolA,
			                  A.rowptr(),
			                  A.colind(),
			                  A.values(),
			                  B.toDense(),
			                  yin,
			                  xout);
		} else {
			// B is sparse
			csr_kron_mult(transA,
			              transB,
			              nrowA,
			              ncolA,
			              A.rowptr(),
		                  A.colind(),
		                  A.values(),
			              nrowB,
			              ncolB,
			              B.rowptr(),
		                  B.colind(),
		                  B.values(),
			              yin,
			              xout);
		};
	};
} // kron_mult

} // namespace Dmrg
#endif // MATRIXDENSEORSPARSE_H
