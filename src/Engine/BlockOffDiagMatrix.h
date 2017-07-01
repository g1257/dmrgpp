#ifndef BLOCKOFFDIAGMATRIX_H
#define BLOCKOFFDIAGMATRIX_H
#include "CrsMatrix.h"
#include "BlockDiagonalMatrix.h"

namespace Dmrg {

template<typename MatrixBlockType>
class BlockOffDiagMatrix {

	typedef typename MatrixBlockType::value_type ComplexOrRealType;
	typedef PsimagLite::CrsMatrix<ComplexOrRealType> SparseMatrixType;
	typedef BlockDiagonalMatrix<MatrixBlockType> BlockDiagonalMatrixType;
	typedef typename PsimagLite::Real<ComplexOrRealType>::Type RealType;

public:

	BlockOffDiagMatrix()
	{}

	BlockOffDiagMatrix(const SparseMatrixType&)
	{
		err("BlockOffDiagMatrix: ctor from sparse\n");
	}

	void fromBlockDiagonal(const BlockDiagonalMatrixType&)
	{
		err("BlockOffDiagMatrix: fromBlockDiagonal unimplemented\n");
	}

	SizeType rows() const
	{
		err("BlockOffDiagMatrix: rows unimplemented\n");
		return 0;
	}

	SizeType cols() const
	{
		err("BlockOffDiagMatrix: cols unimplemented\n");
		return 0;
	}

	SizeType blocks() const
	{
		err("BlockOffDiagMatrix: blocks unimplemented\n");
		return 0;
	}

	void toFull(MatrixBlockType&) const
	{
		err("BlockOffDiagMatrix: toFull unimplemented\n");
	}

	void fromFull(const MatrixBlockType&)
	{
		err("BlockOffDiagMatrix: fromFull unimplemented\n");
	}

	BlockOffDiagMatrix& operator*=(ComplexOrRealType val)
	{
		err("BlockOffDiagMatrix: operator*= unimplemented\n");
		return *this;
	}

	friend bool isTheIdentity(const BlockOffDiagMatrix&)
	{
		err("BlockOffDiagMatrix: isTheIdentity unimplemented\n");
		return false;
	}

	friend void outerProduct(BlockOffDiagMatrix&,
	                         const BlockOffDiagMatrix&,
	                         int nout,
	                         const typename PsimagLite::Vector<RealType>::Type& signs,
	                         bool order)
	{
		err("BlockOffDiagMatrix: externalProduct unimplemented\n");
	}

	friend void transposeConjugate(BlockOffDiagMatrix&,
	                               const BlockOffDiagMatrix&)
	{
		err("BlockOffDiagMatrix: transposeConjugate unimplemented\n");
	}

	friend std::istream& operator>>(std::istream& is,
	                                const BlockOffDiagMatrix& b)
	{
		err("BlockOffDiagMatrix: operator>> unimplemented\n");
		return is;
	}

	friend std::ostream& operator<<(std::ostream& os,
	                                const BlockOffDiagMatrix& b)
	{
		err("BlockOffDiagMatrix: operator>> unimplemented\n");
		return os;
	}
};
}
#endif // BLOCKOFFDIAGMATRIX_H
