#ifndef WFTACCELPATCHES_H
#define WFTACCELPATCHES_H
#include "Matrix.h"
#include "BLAS.h"
#include "ProgramGlobals.h"
#include "MatrixVectorKron/GenIjPatch.h"
#include "BlockDiagWf.h"

namespace Dmrg {

template<typename WaveFunctionTransfBaseType>
class WftAccelPatches {

	typedef typename WaveFunctionTransfBaseType::DmrgWaveStructType DmrgWaveStructType;
	typedef typename WaveFunctionTransfBaseType::WftOptionsType WftOptionsType;
	typedef typename WaveFunctionTransfBaseType::VectorWithOffsetType VectorWithOffsetType;
	typedef typename WaveFunctionTransfBaseType::VectorSizeType VectorSizeType;
	typedef typename DmrgWaveStructType::LeftRightSuperType LeftRightSuperType;
	typedef typename VectorWithOffsetType::VectorType VectorType;
	typedef typename VectorType::value_type ComplexOrRealType;
	typedef typename DmrgWaveStructType::BasisWithOperatorsType BasisWithOperatorsType;
	typedef typename BasisWithOperatorsType::SparseMatrixType SparseMatrixType;
	typedef typename WaveFunctionTransfBaseType::PackIndicesType PackIndicesType;
	typedef typename DmrgWaveStructType::BlockDiagonalMatrixType BlockDiagonalMatrixType;
	typedef typename BlockDiagonalMatrixType::BuildingBlockType MatrixType;
	typedef GenIjPatch<LeftRightSuperType> GenIjPatchType;
	typedef BlockDiagWf<GenIjPatchType, VectorWithOffsetType> BlockDiagWfType;

public:

	WftAccelPatches(const DmrgWaveStructType& dmrgWaveStruct,
	                const WftOptionsType& wftOptions)
	    : dmrgWaveStruct_(dmrgWaveStruct), wftOptions_(wftOptions)
	{}

	void operator()(VectorWithOffsetType& psiDest,
	                SizeType iNew,
	                const VectorWithOffsetType& psiSrc,
	                SizeType iOld,
	                const LeftRightSuperType& lrs,
	                const VectorSizeType& nk,
	                typename ProgramGlobals::DirectionEnum dir) const
	{
		char charLeft = (dir == ProgramGlobals::DirectionEnum::EXPAND_SYSTEM) ? 'C' : 'N';
		char charRight = (dir == ProgramGlobals::DirectionEnum::EXPAND_SYSTEM) ? 'T' : 'N';

		BlockDiagWfType psi(psiSrc,
		                    iOld,
		                    dmrgWaveStruct_.lrs());

		psi.transform(charLeft,
		              charRight,
		              dmrgWaveStruct_.getTransform(ProgramGlobals::SysOrEnvEnum::SYSTEM),
		              dmrgWaveStruct_.getTransform(ProgramGlobals::SysOrEnvEnum::ENVIRON),
		              wftOptions_.gemmRnb,
		              wftOptions_.threadsForGemmR);

		psi.toVectorWithOffsets(psiDest, iNew, lrs, nk, dir);
	}

private:

	const DmrgWaveStructType& dmrgWaveStruct_;
	const WftOptionsType& wftOptions_;
};
}
#endif // WFTACCELPATCHES_H
