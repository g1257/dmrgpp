#ifndef WFTACCELWITHTEMP_H
#define WFTACCELWITHTEMP_H
#include "Matrix.h"
#include "BLAS.h"
#include "ProgramGlobals.h"

namespace Dmrg {

template<typename WaveFunctionTransfBaseType, typename MatrixOrIdentityType>
class WftAccelWithTemp {

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

	WftAccelWithTemp(const DmrgWaveStructType& dmrgWaveStruct,
	                 const WftOptionsType& wftOptions)
	    : dmrgWaveStruct_(dmrgWaveStruct), wftOptions_(wftOptions)
	{}

	void systemFromInfinite(VectorType& dest,
	                        SizeType destOffset,
	                        const VectorType& psiV,
	                        SizeType offset,
	                        const LeftRightSuperType& lrs,
	                        const VectorSizeType& nk,
	                        const SparseMatrixType& ws,
	                        const SparseMatrixType& we) const
	{
		SizeType volumeOfNk = DmrgWaveStructType::volumeOf(nk);
		SizeType nip = lrs.left().permutationInverse().size()/volumeOfNk;
		SizeType nalpha = lrs.left().permutationInverse().size();
		SizeType ni = dmrgWaveStruct_.lrs.right().size()/volumeOfNk;
		MatrixOrIdentityType weRef(wftOptions_.twoSiteDmrg && ni>volumeOfNk, we);

		assert(nip==dmrgWaveStruct_.ws.cols());

		MatrixType temporary;
		buildTemporarySystem(temporary, psiV, offset, ws);

		PackIndicesType pack1(nalpha);
		PackIndicesType pack2(nip);
		for (SizeType x=0;x<dest.size();x++) {
			SizeType isn,jen;
			pack1.unpack(isn,jen,(SizeType)lrs.super().permutation(x+destOffset));
			SizeType is = 0;
			SizeType jpl = 0;
			pack2.unpack(is,jpl,(SizeType)lrs.left().permutation(isn));
			ComplexOrRealType sum = 0;
			for (SizeType k2=weRef.getRowPtr(jen);k2<weRef.getRowPtr(jen+1);k2++) {
				int jpr = weRef.getColOrExit(k2);
				// jpr < 0 could be due to an m smaller than h, the Hilbert size of one site
				// this is checked against elsewhere
				assert(jpr >= 0);
				SizeType jp = dmrgWaveStruct_.lrs.right().permutationInverse(jpl + jpr*volumeOfNk);
				sum += temporary(is, jp)*weRef.getValue(k2);
			}

			dest[x] += sum;
		}
	}

	void environFromInfinite(VectorWithOffsetType& psiDest,
	                         const VectorWithOffsetType& psiSrc,
	                         const LeftRightSuperType& lrs,
	                         SizeType i0,
	                         const VectorSizeType& nk) const
	{
		SizeType volumeOfNk = DmrgWaveStructType::volumeOf(nk);
		SizeType nip = lrs.super().permutationInverse().size()/
		        lrs.right().permutationInverse().size();

		SizeType nip2 = dmrgWaveStruct_.lrs.left().size()/volumeOfNk;

		assert(lrs.left().permutationInverse().size()==volumeOfNk ||
		       lrs.left().permutationInverse().size()==dmrgWaveStruct_.ws.rows());
		assert(lrs.right().permutationInverse().size()/volumeOfNk==dmrgWaveStruct_.we.cols());

		SizeType start = psiDest.offset(i0);
		SizeType total = psiDest.effectiveSize(i0);

		SparseMatrixType ws;
		dmrgWaveStruct_.ws.toSparse(ws);
		SparseMatrixType we;
		dmrgWaveStruct_.we.toSparse(we);
		SparseMatrixType weT;
		transposeConjugate(weT,we);

		MatrixOrIdentityType wsRef2(wftOptions_.twoSiteDmrg && nip>volumeOfNk,ws);

		MatrixType temporary;
		buildTemporaryEnviron(temporary, psiSrc, i0, we);

		PackIndicesType pack1(nip);
		PackIndicesType pack2(volumeOfNk);
		for (SizeType x=0;x<total;x++) {
			SizeType ip = 0;
			SizeType beta = 0;
			pack1.unpack(ip,beta,(SizeType)lrs.super().permutation(x+start));
			SizeType kp = 0;
			SizeType jp = 0;
			pack2.unpack(kp, jp, lrs.right().permutation(beta));
			for (SizeType k3=wsRef2.getRowPtr(ip);k3<wsRef2.getRowPtr(ip+1);k3++) {
				int ip2 = wsRef2.getColOrExit(k3);
				if (ip2<0) continue;
				SizeType alpha = dmrgWaveStruct_.lrs.left().permutationInverse(ip2+kp*nip2);
				psiDest.fastAccess(i0,x) += temporary(jp, alpha)*wsRef2.getValue(k3);
			}
		}
	}

private:

	void buildTemporarySystem(MatrixType& temporary,
	                          const VectorType& psiV,
	                          SizeType offset,
	                          const SparseMatrixType& ws) const
	{
		SizeType nalpha=dmrgWaveStruct_.lrs.left().permutationInverse().size();
		temporary.resize(nalpha, dmrgWaveStruct_.lrs.super().size()/nalpha);
		temporary.setTo(0.0);

		PackIndicesType packSuper(nalpha);
		SizeType total = psiV.size();

		for (SizeType y = 0; y < total; ++y) {
			SizeType ip = 0;
			SizeType jp = 0;
			packSuper.unpack(ip, jp, dmrgWaveStruct_.lrs.super().permutation(y + offset));
			SizeType start = ws.getRowPtr(ip);
			SizeType end = ws.getRowPtr(ip + 1);
			for (SizeType k = start;k < end; k++) {
				SizeType is = ws.getCol(k);
				temporary(is, jp) += ws.getValue(k)*psiV[y];
			}
		}
	}

	void buildTemporaryEnviron(MatrixType& temporary,
	                           const VectorWithOffsetType& psiV,
	                           SizeType i0,
	                           const SparseMatrixType& we) const
	{
		SizeType ni=dmrgWaveStruct_.lrs.left().size();
		temporary.resize(we.cols(), ni);
		temporary.setTo(0.0);

		PackIndicesType packSuper(ni);
		SizeType total = psiV.effectiveSize(i0);
		SizeType offset = psiV.offset(i0);

		for (SizeType x = 0; x < total; ++x) {
			SizeType alpha = 0;
			SizeType jp2 = 0;
			packSuper.unpack(alpha, jp2, dmrgWaveStruct_.lrs.super().permutation(x + offset));
			SizeType start = we.getRowPtr(jp2);
			SizeType end = we.getRowPtr(jp2+ 1);
			for (SizeType k = start; k < end; k++) {
				SizeType jp = we.getCol(k);
				temporary(jp, alpha) += psiV.fastAccess(i0, x)*we.getValue(k);
			}
		}
	}

	const DmrgWaveStructType& dmrgWaveStruct_;
	const WftOptionsType& wftOptions_;
};
}
#endif // WFTACCELWITHTEMP_H
