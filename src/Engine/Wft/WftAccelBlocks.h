#ifndef WFTACCELBLOCKS_H
#define WFTACCELBLOCKS_H
#include "Matrix.h"
#include "BLAS.h"

namespace Dmrg {

template<typename WaveFunctionTransfBaseType>
class WftAccelBlocks {

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

	WftAccelBlocks(const DmrgWaveStructType& dmrgWaveStruct,
	               const WftOptionsType& wftOptions)
	    : dmrgWaveStruct_(dmrgWaveStruct), wftOptions_(wftOptions)
	{}

	void environFromInfinite(VectorWithOffsetType& psiDest,
	                         SizeType i0,
	                         const VectorWithOffsetType& psiSrc,
	                         SizeType i0src,
	                         const LeftRightSuperType& lrs,
	                         const VectorSizeType& nk) const
	{
		SizeType volumeOfNk = DmrgWaveStructType::volumeOf(nk);
		SizeType jp2size = dmrgWaveStruct_.lrs.super().size();
		SizeType i2psize = dmrgWaveStruct_.lrs.left().permutationInverse().size()/volumeOfNk;

		VectorMatrixType psi(volumeOfNk);
		for (SizeType kp = 0; kp < volumeOfNk; ++kp) {
			psi[kp].resize(i2psize, jp2size);
			psi[kp].setTo(0.0);
		}

		environPreparePsi(psi, psiSrc, i0src, volumeOfNk);
		MatrixType ws;
		dmrgWaveStruct_.ws.toDense(ws);

		MatrixType we;
		dmrgWaveStruct_.we.toDense(we);

		SizeType jpsize = we.cols();
		SizeType ipsize = ws.rows();
		MatrixType tmp(i2psize, jpsize);
		VectorMatrixType result(volumeOfNk);

		for (SizeType kp = 0; kp < volumeOfNk; ++kp) {
			result[kp].resize(ipsize, jpsize);
			result[kp].setTo(0.0);
			SizeType ip = result[kp].rows();
			tmp.setTo(0.0);

			psimag::BLAS::GEMM('N',
			                   'N',
			                   i2psize,
			                   jpsize,
			                   jp2size,
			                   1.0,
			                   &((psi[kp])(0,0)),
			                   i2psize,
			                   &(we(0,0)),
			                   jp2size,
			                   0.0,
			                   &(tmp(0,0)),
			                   i2psize);

			psimag::BLAS::GEMM('N',
			                   'N',
			                   ipsize,
			                   jpsize,
			                   i2psize,
			                   1.0,
			                   &(ws(0,0)),
			                   ip,
			                   &(tmp(0,0)),
			                   i2psize,
			                   0.0,
			                   &((result[kp])(0,0)),
			                   ipsize);
		}

		environCopyOut(psiDest, i0, result, lrs, volumeOfNk);
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

private:

	void environPreparePsi(VectorMatrixType& psi,
	                       const VectorWithOffsetType& psiSrc,
	                       SizeType i0src,
	                       SizeType volumeOfNk) const
	{
		SizeType total = psiSrc.effectiveSize(i0src);
		SizeType offset = psiSrc.offset(i0src);
		PackIndicesType packSuper(dmrgWaveStruct_.lrs.left().size());
		PackIndicesType packLeft(dmrgWaveStruct_.lrs.left().size()/volumeOfNk);

		for (SizeType x = 0; x < total; ++x) {
			SizeType alpha = 0;
			SizeType jp2 = 0;
			packSuper.unpack(alpha, jp2, dmrgWaveStruct_.lrs.super().permutation(x + offset));
			SizeType ip2 = 0;
			SizeType kp = 0;
			packLeft.unpack(ip2, kp, dmrgWaveStruct_.lrs.left().permutation(alpha));
			psi[kp](ip2, jp2) += psiSrc.fastAccess(i0src, x);
		}
	}

	void environCopyOut(VectorWithOffsetType& psiDest,
	                    SizeType i0,
	                    const VectorMatrixType& result,
	                    const LeftRightSuperType& lrs,
	                    SizeType volumeOfNk) const
	{
		SizeType nip = lrs.super().permutationInverse().size()/
		        lrs.right().permutationInverse().size();
		PackIndicesType pack1(nip);
		PackIndicesType pack2(volumeOfNk);
		SizeType total = psiDest.effectiveSize(i0);
		SizeType start = psiDest.offset(i0);

		for (SizeType x = 0; x < total; ++x) {
			SizeType ip = 0;
			SizeType beta = 0;
			pack1.unpack(ip,beta,(SizeType)lrs.super().permutation(x+start));
			SizeType kp = 0;
			SizeType jp = 0;
			pack2.unpack(kp,jp,(SizeType)lrs.right().permutation(beta));
			psiDest.fastAccess(i0, x) += result[kp](ip, jp);
		}
	}

	const DmrgWaveStructType& dmrgWaveStruct_;
	const WftOptionsType& wftOptions_;
};
}
#endif // WFTACCELBLOCKS_H
