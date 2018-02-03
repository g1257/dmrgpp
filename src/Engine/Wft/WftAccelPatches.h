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
	typedef typename WaveFunctionTransfBaseType::WftOptions WftOptionsType;
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
	                const VectorWithOffsetType& psiSrc,
	                const LeftRightSuperType& lrs,
	                SizeType ii,
	                const VectorSizeType& nk,
	                typename ProgramGlobals::DirectionEnum dir) const
	{
		char charLeft = (dir == ProgramGlobals::EXPAND_SYSTEM) ? 'C' : 'N';
		char charRight = charLeft;
		const BlockDiagonalMatrixType& transformLeft = dmrgWaveStruct_.ws;
		const BlockDiagonalMatrixType& transformRight = dmrgWaveStruct_.we;
		BlockDiagWfType psi(psiSrc,
		                         psiSrc.sector(ii),
		                         dmrgWaveStruct_.lrs);

		psi.transform(charLeft, charRight, transformLeft, transformRight);
		psi.toVectorWithOffsets(psiDest);
	}

private:

	const DmrgWaveStructType& dmrgWaveStruct_;
	const WftOptionsType& wftOptions_;
};
}
#endif // WFTACCELPATCHES_H
