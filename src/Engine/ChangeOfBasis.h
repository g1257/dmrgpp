#ifndef DMRG_CHANGEOFBASIS_H
#define DMRG_CHANGEOFBASIS_H
#include "BlockDiagonalMatrix.h"
#include "BlockOffDiagMatrix.h"

namespace Dmrg {

template<typename SparseMatrixType, typename MatrixType>
class ChangeOfBasis {

public:

	typedef BlockDiagonalMatrix<MatrixType> BlockDiagonalMatrixType;
	typedef BlockOffDiagMatrix<MatrixType> BlockOffDiagMatrixType;

	void update(const BlockDiagonalMatrixType& transform)
	{
		transform_ = transform;
	}

	void operator()(SparseMatrixType &v) const
	{
		BlockOffDiagMatrixType vBlocked(v, transform_.offsetsRows());
		vBlocked.transform(transform_);
		vBlocked.toSparse(v);
	}

	static void changeBasis(SparseMatrixType &v,
	                        const BlockDiagonalMatrixType& ftransform1)
	{
		BlockOffDiagMatrixType vBlocked(v, ftransform1.offsetsRows());
		vBlocked.transform(ftransform1);
		vBlocked.toSparse(v);
	}

private:

	BlockDiagonalMatrixType transform_;
}; // class ChangeOfBasis
} // namespace Dmrg
#endif // DMRG_CHANGEOFBASIS_H
