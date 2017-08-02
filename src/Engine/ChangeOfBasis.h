#ifndef DMRG_CHANGEOFBASIS_H
#define DMRG_CHANGEOFBASIS_H
#include "BlockDiagonalMatrix.h"

namespace Dmrg {

template<typename SparseMatrixType, typename MatrixType>
class ChangeOfBasis {

public:

	typedef BlockDiagonalMatrix<MatrixType> BlockDiagonalMatrixType;

	void update(const BlockDiagonalMatrixType& transform)
	{
		transform.toSparse(transform_);
		transposeConjugate(transformT_,transform_);
	}

	void update(const MatrixType& v2)
	{
		fullMatrixToCrsMatrix(transform_, v2);
		transposeConjugate(transformT_,transform_);
	}

	void operator()(SparseMatrixType &v) const
	{
		SparseMatrixType tmp;
		multiply(tmp,v,transform_);
		multiply(v,transformT_,tmp);
	}

	static void changeBasis(SparseMatrixType &v,
	                        const BlockDiagonalMatrixType& ftransform1)
	{
		SparseMatrixType ftransform;
		ftransform1.toSparse(ftransform);
		SparseMatrixType ftransformT;
		transposeConjugate(ftransformT,ftransform);
		SparseMatrixType tmp;
		multiply(tmp,v,ftransform);
		multiply(v,ftransformT,tmp);
	}

private:

	SparseMatrixType transform_;
	SparseMatrixType transformT_;
}; // class ChangeOfBasis
} // namespace Dmrg
#endif // DMRG_CHANGEOFBASIS_H
