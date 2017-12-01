#ifndef WFTACCELBLOCKS_H
#define WFTACCELBLOCKS_H

namespace Dmrg {

template<typename WaveFunctionTransfBaseType>
class WftAccelBlocks {

	typedef typename WaveFunctionTransfBaseType::DmrgWaveStructType DmrgWaveStructType;
	typedef typename WaveFunctionTransfBaseType::VectorWithOffsetType VectorWithOffsetType;
	typedef typename WaveFunctionTransfBaseType::VectorSizeType VectorSizeType;
	typedef typename DmrgWaveStructType::LeftRightSuperType LeftRightSuperType;
	typedef typename VectorWithOffsetType::VectorType VectorType;
	typedef typename DmrgWaveStructType::BasisWithOperatorsType BasisWithOperatorsType;
	typedef typename BasisWithOperatorsType::SparseMatrixType SparseMatrixType;

public:

	void environFromInfinite(VectorWithOffsetType& psiDest,
	                         const VectorWithOffsetType& psiSrc,
	                         const LeftRightSuperType& lrs,
	                         SizeType i0,
	                         const VectorSizeType& nk) const
	{
//		for (SizeType k = 0; k < kmax; ++k) {
//			blocksPreparePsi(psi[k], psiSrc, k);
//			PsimagLite::BLAS::GEMM('N',
//			                       'N',
//			                       i2psize,
//			                       jpsize,
//			                       jp2size,
//			                       1.0,
//			                       psi[k],
//			                       i2psize,
//			                       we,
//			                       j2psize,
//			                       0.0,
//			                       tmp,
//			                       i2psize);

//			PsimagLite::BLAS::GEMM('N',
//			                       'N',
//			                       ipSize,
//			                       jpsize,
//			                       i2psize,
//			                       1.0,
//			                       ws,
//			                       ip,
//			                       tmp,
//			                       i2psize,
//			                       0.0,
//			                       result[k],
//			                       ipsize);
//		}

//		blocksCopyOut(psiDest, i0, result);
		err("environFromInfinite not yet implemented\n");
	}

	void systemFromInfinite(VectorType& dest,
	                        SizeType destOffset,
	                        const VectorType& psiV,
	                        SizeType offset,
	                        const LeftRightSuperType& lrs,
	                        const VectorSizeType& nk,
	                        const SparseMatrixType& wsT,
	                        const SparseMatrixType& we) const
	{
		err("systemFromInfinite not yet implemented\n");
	}

};
}
#endif // WFTACCELBLOCKS_H
