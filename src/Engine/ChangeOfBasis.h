#ifndef DMRG_CHANGEOFBASIS_H
#define DMRG_CHANGEOFBASIS_H
#include "BlockDiagonalMatrix.h"
#include "BlockOffDiagMatrix.h"

namespace Dmrg {

template<typename SparseMatrixType, typename MatrixType>
class ChangeOfBasis {

#ifdef JZ_SYMMETRY
	static const bool jzSymmetry_ = true;
#else
	static const bool jzSymmetry_ = false;
#endif

public:

	typedef BlockDiagonalMatrix<MatrixType> BlockDiagonalMatrixType;
	typedef BlockOffDiagMatrix<MatrixType> BlockOffDiagMatrixType;

	ChangeOfBasis()
	{
		if (!jzSymmetry_) return;
		std::cout<<"JZ_SYMMETRY in use\n";
	}

	void update(const BlockDiagonalMatrixType& transform)
	{
		if (!jzSymmetry_) {
			transform_ = transform;
			return;
		}

		assert(jzSymmetry_);
		transform.toSparse(oldT_);
		transposeConjugate(oldTtranspose_, oldT_);
	}

	void operator()(SparseMatrixType &v) const
	{
		if (!jzSymmetry_) {
			BlockOffDiagMatrixType vBlocked(v, transform_.offsetsRows());
			vBlocked.transform(transform_);
			vBlocked.toSparse(v);
			return;
		}

		assert(jzSymmetry_);
		SparseMatrixType tmp;
		multiply(tmp, v, oldT_);
		multiply(v, oldTtranspose_, tmp);
	}

	static void changeBasis(SparseMatrixType &v,
	                        const BlockDiagonalMatrixType& ftransform1)
	{
		if (!jzSymmetry_) {
			BlockOffDiagMatrixType vBlocked(v, ftransform1.offsetsRows());
			vBlocked.transform(ftransform1);
			vBlocked.toSparse(v);
			return;
		}

		assert(jzSymmetry_);
		SparseMatrixType ftransform;
		ftransform1.toSparse(ftransform);
		SparseMatrixType ftransformT;
		transposeConjugate(ftransformT,ftransform);
		SparseMatrixType tmp;
		multiply(tmp,v,ftransform);
		multiply(v,ftransformT,tmp);
	}

private:

	BlockDiagonalMatrixType transform_;
	SparseMatrixType oldT_;
	SparseMatrixType oldTtranspose_;
}; // class ChangeOfBasis
} // namespace Dmrg
#endif // DMRG_CHANGEOFBASIS_H
