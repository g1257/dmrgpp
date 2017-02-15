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
	    : isDense_(sparse.nonZero() > static_cast<int>(0.5*sparse.row()*sparse.col())),
	      rows_(sparse.row()),
	      cols_(sparse.col())
	{
		if (!isDense_) {
			rowptr_.resize(rows_ + 1, 0);
			for (SizeType k = 0; k < rows_ + 1; ++k)
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
			values_.resize(rows_*cols_, 0.0);
			for (SizeType i = 0; i < rows_; ++i) {
				for (int k = sparse.getRowPtr(i); k < sparse.getRowPtr(i+1); ++k)
					values_[i + sparse.getCol(k)*rows_] = sparse.getValue(k);
			}
		}
	}

	bool isDense() const { return isDense_; }

	SizeType rows() const { return rows_; }

	SizeType cols() const { return cols_; }

	const VectorType& values() const { return values_; }

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
			PsimagLite::Matrix<ComplexOrRealType> m2(rows_, cols_);
			for (SizeType i = 0; i < rows_; ++i)
				for (SizeType j = 0; j < cols_; ++j)
					m2(i,j) = values_[i + j*rows_];

			return SparseMatrixType(m2);
		}

		SparseMatrixType m;
		m.resize(rows_, cols_, colind_.size());
		SizeType counter = 0;
		for (SizeType i = 0; i < rows_; ++i) {
			m.setRow(i, counter);
			for (int k = rowptr_[i]; k < rowptr_[i+1]; ++k)
				m.setCol(counter++,colind_[k]);
		}

		m.setRow(rows_, counter);
		m.checkValidity();
		return m;
	}

private:

	bool isDense_;
	SizeType rows_;
	SizeType cols_;
	VectorIntType rowptr_;
	VectorIntType colind_;
	VectorType values_;
}; // class MatrixDenseOrSparse

template<typename SparseMatrixType>
void kronMult(typename PsimagLite::Vector<typename SparseMatrixType::value_type>::Type& xout,
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

	const typename SparseMatrixType::value_type *aval = &(A.values()[0]);
	const typename SparseMatrixType::value_type *bval = &(B.values()[0]);

	const bool isDenseA = A.isDense();
	const bool isDenseB = B.isDense();

	if (isDenseA) {
		if (isDenseB) {
			den_kron_mult(transA,
			              transB,
			              nrowA,
			              ncolA,
			              aval,
			              nrowB,
			              ncolB,
			              bval,
			              yin,
			        &(xout[0]));
		} else  {
			// B is sparse
			const int *browptr = &(B.rowptr()[0]);
			const int *bcol = &(B.colind()[0]);

			den_csr_kron_mult(transA,
			                  transB,
			                  nrowA,
			                  ncolA,
			                  aval,
			                  nrowB,
			                  ncolB,
			                  browptr,
			                  bcol,
			                  bval,
			                  yin,
			        &(xout[0]));
		}
	} else {
		// A is sparse
		const int *arowptr = &(A.rowptr()[0]);
		const int *acol = &(A.colind()[0]);

		if (isDenseB) {
			csr_den_kron_mult(transA,
			                  transB,
			                  nrowA,
			                  ncolA,
			                  arowptr,
			                  acol,
			                  aval,
			                  nrowB,
			                  ncolB,
			                  bval,
			                  yin,
			        &(xout[0]));
		} else {
			// B is sparse
			const int *browptr = &(B.rowptr()[0]);
			const int *bcol = &(B.colind()[0]);

			csr_kron_mult(transA,
			              transB,
			              nrowA,
			              ncolA,
			              arowptr,
			              acol,
			              aval,
			              nrowB,
			              ncolB,
			              browptr,
			              bcol,
			              bval,
			              yin,
			        &(xout[0]));
		};
	};
} // kron_mult

} // namespace Dmrg
#endif // MATRIXDENSEORSPARSE_H
