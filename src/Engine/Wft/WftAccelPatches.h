#ifndef WFTACCELPATCHES_H
#define WFTACCELPATCHES_H
#include "Matrix.h"
#include "BLAS.h"
#include "ProgramGlobals.h"

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
	typedef PsimagLite::Matrix<ComplexOrRealType> MatrixType;
	typedef typename PsimagLite::Vector<MatrixType>::Type VectorMatrixType;
	typedef typename WaveFunctionTransfBaseType::PackIndicesType PackIndicesType;

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
		SizeType qn = psiDest.qn(ii);
		SizeType iOld = findIold(psiSrc, qn);
		SizeType i0 = psiDest.sector(ii);
		VectorType psiDestOneSector;
		psiDest.extract(psiDestOneSector, i0);

		VectorType psiSrcOneSector;
		psiSrc.extract(psiSrcOneSector, iOld);

		psiDest.setDataInSector(psiDestOneSector, i0);
	}

	static SizeType findIold(const VectorWithOffsetType& psiSrc,
	                  SizeType qn)
	{
		SizeType sectors = psiSrc.sectors();
		for (SizeType i = 0; i < sectors; ++i)
			if (psiSrc.qn(i) == qn)
				return psiSrc.sector(i);

		err("WaveFunctionTransfLocal::findIold(): Cannot find sector in old vector\n");
		throw PsimagLite::RuntimeError("UNREACHABLE\n");
	}

private:

	const DmrgWaveStructType& dmrgWaveStruct_;
	const WftOptionsType& wftOptions_;
};
}
#endif // WFTACCELPATCHES_H
