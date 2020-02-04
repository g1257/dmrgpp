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

	void operator()(OperatorStorageType& v,
	                SizeType gemmRnb,
	                SizeType threadsForGemmR) const
	{
		if (!ProgramGlobals::oldChangeOfBasis) {
			BlockOffDiagMatrixType vBlocked(v.getCRS(), transform_.offsetsRows());
			vBlocked.transform(transform_, gemmRnb, threadsForGemmR);
			vBlocked.toSparse(v.getCRSNonConst());
			return;
		}

		assert(ProgramGlobals::oldChangeOfBasis);
		v.rotate(oldTtranspose_, oldT_);
	}

	static void changeBasis(OperatorStorageType &v,
	                        const BlockDiagonalMatrixType& ftransform1,
	                        SizeType gemmRnb,
	                        SizeType threadsForGemmR)
	{
		if (!v.justCRS())
			err("changeBasis: operatorstorage not justCRS\n");

		if (!ProgramGlobals::oldChangeOfBasis) {
			BlockOffDiagMatrixType vBlocked(v.getCRS(), ftransform1.offsetsRows());
			vBlocked.transform(ftransform1, gemmRnb, threadsForGemmR);
			vBlocked.toSparse(v.getCRSNonConst());
			return;
		}

		assert(ProgramGlobals::oldChangeOfBasis);
		SparseMatrixType ftransform;
		ftransform1.toSparse(ftransform);
		SparseMatrixType ftransformT;
		transposeConjugate(ftransformT,ftransform);
		v.rotate(ftransformT, ftransform);
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
