#ifndef DMRG_CHANGEOFBASIS_H
#define DMRG_CHANGEOFBASIS_H
#include "BlockDiagonalMatrix.h"
#include "BlockOffDiagMatrix.h"
#include "ProgramGlobals.h"

namespace Dmrg {

template<typename OperatorStorageType, typename MatrixType>
class ChangeOfBasis {

public:

	typedef BlockDiagonalMatrix<MatrixType> BlockDiagonalMatrixType;
	typedef BlockOffDiagMatrix<MatrixType> BlockOffDiagMatrixType;
	typedef typename OperatorStorageType::value_type ComplexOrRealType;
	typedef PsimagLite::CrsMatrix<ComplexOrRealType> SparseMatrixType;

	ChangeOfBasis()
	{
		if (!ProgramGlobals::oldChangeOfBasis) return;
		std::cout<<"Old ChangeOfBasis in use\n";
	}

	void update(const BlockDiagonalMatrixType& transform)
	{
		if (!ProgramGlobals::oldChangeOfBasis) {
			transform_ = transform;
			return;
		}

		assert(ProgramGlobals::oldChangeOfBasis);
		transform.toSparse(oldT_);
		transposeConjugate(oldTtranspose_, oldT_);
	}

#ifndef USE_OPERATOR_STORAGE
	void operator()(SparseMatrixType &v) const
	{
		if (!ProgramGlobals::oldChangeOfBasis) {
			BlockOffDiagMatrixType vBlocked(v, transform_.offsetsRows());
			vBlocked.transform(transform_);
			vBlocked.toSparse(v);
			return;
		}

		assert(ProgramGlobals::oldChangeOfBasis);
		SparseMatrixType tmp;
		multiply(tmp, v, oldT_);
		multiply(v, oldTtranspose_, tmp);
	}
#else
	void operator()(OperatorStorageType &v) const
	{
		err("Needs to be written for OperatorStorageType\n");
	}
#endif

	static void changeBasis(SparseMatrixType &v,
	                        const BlockDiagonalMatrixType& ftransform1)
	{
		if (!ProgramGlobals::oldChangeOfBasis) {
			BlockOffDiagMatrixType vBlocked(v, ftransform1.offsetsRows());
			vBlocked.transform(ftransform1);
			vBlocked.toSparse(v);
			return;
		}

		assert(ProgramGlobals::oldChangeOfBasis);
		SparseMatrixType ftransform;
		ftransform1.toSparse(ftransform);
		SparseMatrixType ftransformT;
		transposeConjugate(ftransformT,ftransform);
		SparseMatrixType tmp;
		multiply(tmp,v,ftransform);
		multiply(v,ftransformT,tmp);
	}

	void clear()
	{
		transform_.clear();
		oldT_.clear();
		oldTtranspose_.clear();
	}

private:

	BlockDiagonalMatrixType transform_;
	SparseMatrixType oldT_;
	SparseMatrixType oldTtranspose_;
}; // class ChangeOfBasis
} // namespace Dmrg
#endif // DMRG_CHANGEOFBASIS_H
