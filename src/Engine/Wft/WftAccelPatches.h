#ifndef WFTACCELPATCHES_H
#define WFTACCELPATCHES_H
#include "BlockDiagWf.h"
#include "Matrix.h"
#include "MatrixVectorKron/GenIjPatch.h"
#include "ProgramGlobals.h"

namespace Dmrg {

template <typename WaveFunctionTransfBaseType> class WftAccelPatches {

	using DmrgWaveStructType      = typename WaveFunctionTransfBaseType::DmrgWaveStructType;
	using WftOptionsType          = typename WaveFunctionTransfBaseType::WftOptionsType;
	using VectorWithOffsetType    = typename WaveFunctionTransfBaseType::VectorWithOffsetType;
	using VectorSizeType          = typename WaveFunctionTransfBaseType::VectorSizeType;
	using OneSiteSpacesType       = typename WaveFunctionTransfBaseType::OneSiteSpacesType;
	using LeftRightSuperType      = typename DmrgWaveStructType::LeftRightSuperType;
	using VectorType              = typename VectorWithOffsetType::VectorType;
	using ComplexOrRealType       = typename VectorType::value_type;
	using BasisWithOperatorsType  = typename DmrgWaveStructType::BasisWithOperatorsType;
	using SparseMatrixType        = typename BasisWithOperatorsType::SparseMatrixType;
	using PackIndicesType         = typename WaveFunctionTransfBaseType::PackIndicesType;
	using BlockDiagonalMatrixType = typename DmrgWaveStructType::BlockDiagonalMatrixType;
	using MatrixType              = typename BlockDiagonalMatrixType::BuildingBlockType;
	using GenIjPatchType          = GenIjPatch<LeftRightSuperType>;
	using BlockDiagWfType
	    = BlockDiagWf<GenIjPatchType, VectorWithOffsetType, OneSiteSpacesType>;

public:

	WftAccelPatches(const DmrgWaveStructType& dmrgWaveStruct, const WftOptionsType& wftOptions)
	    : dmrgWaveStruct_(dmrgWaveStruct)
	    , wftOptions_(wftOptions)
	{ }

	void operator()(VectorWithOffsetType&       psiDest,
	                SizeType                    iNew,
	                const VectorWithOffsetType& psiSrc,
	                SizeType                    iOld,
	                const LeftRightSuperType&   lrs,
	                const OneSiteSpacesType&    oneSiteSpaces) const
	{
		ProgramGlobals::DirectionEnum dir = oneSiteSpaces.direction();
		char charLeft  = (dir == ProgramGlobals::DirectionEnum::EXPAND_SYSTEM) ? 'C' : 'N';
		char charRight = (dir == ProgramGlobals::DirectionEnum::EXPAND_SYSTEM) ? 'T' : 'N';

		BlockDiagWfType psi(psiSrc, iOld, dmrgWaveStruct_.lrs());

		psi.transform(charLeft,
		              charRight,
		              dmrgWaveStruct_.getTransform(ProgramGlobals::SysOrEnvEnum::SYSTEM),
		              dmrgWaveStruct_.getTransform(ProgramGlobals::SysOrEnvEnum::ENVIRON),
		              wftOptions_.gemmRnb,
		              wftOptions_.threadsForGemmR);

		psi.toVectorWithOffsets(psiDest, iNew, lrs, oneSiteSpaces);
	}

private:

	const DmrgWaveStructType& dmrgWaveStruct_;
	const WftOptionsType&     wftOptions_;
};
}
#endif // WFTACCELPATCHES_H
