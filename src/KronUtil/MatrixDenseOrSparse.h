#ifndef MATRIXDENSEORSPARSE_H
#define MATRIXDENSEORSPARSE_H
#include "Vector.h"
#include "KronUtilWrapper.h"
#include "Matrix.h"

namespace Dmrg {

template<typename SparseMatrixType>
class MatrixDenseOrSparse {

public:

	typedef typename SparseMatrixType::value_type value_type;
	typedef typename SparseMatrixType::value_type ComplexOrRealType;
	typedef typename PsimagLite::Real<ComplexOrRealType>::Type RealType;
	typedef typename PsimagLite::Vector<ComplexOrRealType>::Type VectorType;
	typedef PsimagLite::Vector<int>::Type VectorIntType;

	explicit MatrixDenseOrSparse(const SparseMatrixType& sparse,
	                             const RealType& threshold)
	    : isDense_(sparse.nonZeros() > static_cast<SizeType>(threshold*
	                                                         sparse.rows()*
	                                                         sparse.cols())),
	      sparseMatrix_(sparse)
	{
		sparseMatrix_.checkValidity();

		if (isDense_) // A(i,j) at  val[ (i) + (j)*nrow ]
			crsMatrixToFullMatrix(denseMatrix_, sparse);
	}

	explicit MatrixDenseOrSparse(const SizeType nrows,
	                             const SizeType ncols,
	                             bool  isDense_in )
	    : isDense_( isDense_in ), sparseMatrix_(nrows,ncols), denseMatrix_(0,0)
	{
		if (isDense_) {
			denseMatrix_.clear();
			denseMatrix_.resize(nrows,ncols );
		}
		else  {
			sparseMatrix_.resize(nrows,ncols);
		};
	}


	void conjugate()
	{
		SparseMatrixType& nonconst = const_cast<SparseMatrixType&>(sparseMatrix_);
		nonconst.conjugate();
		if (isDense_)
			denseMatrix_.conjugate();
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
		if (isDense_)
			throw PsimagLite::RuntimeError("FATAL: Matrix isn't sparse\n");

		sparseMatrix_.checkValidity();
		return sparseMatrix_;
	}

	bool isZero() const
	{
		return (isDense_) ? denseMatrix_.isZero()  :
		                    sparseMatrix_.isZero();
	}

	SparseMatrixType toSparse() const
	{
		return (isDense_) ? SparseMatrixType(denseMatrix_) : sparse();
	}

	const PsimagLite::CrsMatrix<ComplexOrRealType>& getSparse() const
	{
		assert( !isDense_ );
		return( sparseMatrix_ );
	}

	PsimagLite::CrsMatrix<ComplexOrRealType>& getSparse()
	{
		assert( !isDense_ );
		return( sparseMatrix_ );
	}

	const PsimagLite::Matrix<ComplexOrRealType>& getDense() const
	{
		assert( isDense_ );
		return( denseMatrix_ );
	}

	PsimagLite::Matrix<ComplexOrRealType>& getDense()
	{
		assert( isDense_ );
		return( denseMatrix_ );
	}


private:

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
              const MatrixDenseOrSparse<SparseMatrixType>& B,
              const typename PsimagLite::Real<typename SparseMatrixType::value_type>::Type
              denseFlopDiscount)
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
			              offsetX,
			              denseFlopDiscount);
		} else  {
			// B is sparse
			den_csr_kron_mult(transA,
			                  transB,
			                  A.dense(),
			                  B.sparse(),
			                  yin,
			                  offsetY,
			                  xout,
			                  offsetX,
			                  denseFlopDiscount);
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
			                  offsetX,
			                  denseFlopDiscount);
		} else {
			// B is sparse
			csr_kron_mult(transA,
			              transB,
			              A.sparse(),
			              B.sparse(),
			              yin,
			              offsetY,
			              xout,
			              offsetX,
			              denseFlopDiscount);
		};
	};
} // kron_mult

} // namespace Dmrg
#endif // MATRIXDENSEORSPARSE_H
