#ifndef WFTSPARSETWOSITE_H
#define WFTSPARSETWOSITE_H
#include "Matrix.h"
#include "BLAS.h"
#include "ProgramGlobals.h"

namespace Dmrg {

template<typename WaveFunctionTransfBaseType, typename MatrixOrIdentityType>
class WftSparseTwoSite {

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

	WftSparseTwoSite(VectorWithOffsetType& dest,
	                 SizeType i0,
	                 const VectorWithOffsetType& src,
	                 SizeType iOld,
	                 const DmrgWaveStructType& dmrgWaveStruct,
	                 const WftOptionsType& wftOptions,
	                 const LeftRightSuperType& lrs,
	                 const VectorSizeType& nk,
	                 const SparseMatrixType& wsT,
	                 const SparseMatrixType& we,
	                 SizeType sysOrEnv)
	    : dest_(dest),
	      i0_(i0),
	      src_(src),
	      iOld_(iOld),
	      dmrgWaveStruct_(dmrgWaveStruct),
	      wftOptions_(wftOptions),
	      lrs_(lrs),
	      wsT_(wsT),
	      we_(we),
	      volumeOfNk_(DmrgWaveStructType::volumeOf(nk)),
	      pack1_((sysOrEnv == ProgramGlobals::SYSTEM) ? lrs.left().permutationInverse().size() :
	                                                    lrs.super().permutationInverse().size()/
	                                                    lrs.right().permutationInverse().size()),
	      pack2_((sysOrEnv == ProgramGlobals::SYSTEM) ?  lrs.left().permutationInverse().size()/
	                                                     volumeOfNk_ : volumeOfNk_),
	      sysOrEnv_(sysOrEnv)
	{
		assert(sysOrEnv == ProgramGlobals::SYSTEM || sysOrEnv == ProgramGlobals::ENVIRON);

		if (sysOrEnv == ProgramGlobals::SYSTEM) {
			assert(lrs.left().permutationInverse().size()/volumeOfNk_ ==
			       dmrgWaveStruct_.ws.cols());

		} else {
			assert(lrs.left().permutationInverse().size()==volumeOfNk_ ||
			       lrs.left().permutationInverse().size()==dmrgWaveStruct_.ws.rows());
			assert(lrs.right().permutationInverse().size()/volumeOfNk_ ==
			       dmrgWaveStruct_.we.cols());
		}
	}

	SizeType tasks() const { return dest_.effectiveSize(i0_); }

	void doTask(SizeType x, SizeType)
	{
		SizeType destOffset = dest_.offset(i0_);
		if (sysOrEnv_ == ProgramGlobals::SYSTEM) {
			SizeType isn = 0;
			SizeType jen = 0;
			pack1_.unpack(isn, jen, lrs_.super().permutation(x + destOffset));
			SizeType is = 0;
			SizeType jpl = 0;
			pack2_.unpack(is, jpl, lrs_.left().permutation(isn));

			dest_.fastAccess(i0_, x) += createAux2bFromInfinite(is, jpl, jen);
		} else {
			SizeType ip = 0;
			SizeType beta = 0;
			pack1_.unpack(ip, beta, lrs_.super().permutation(x + destOffset));
			SizeType kp = 0;
			SizeType jp = 0;
			pack2_.unpack(kp, jp, lrs_.right().permutation(beta));
			dest_.fastAccess(i0_, x) += createAux1bFromInfinite(ip, kp, jp);
		}
	}

private:

	ComplexOrRealType createAux2bFromInfinite(SizeType is,
	                                          SizeType jpl,
	                                          SizeType jen) const
	{
		SizeType offset = src_.offset(iOld_);
		SizeType nalpha=dmrgWaveStruct_.lrs.left().permutationInverse().size();
		SizeType ni = dmrgWaveStruct_.lrs.right().size()/volumeOfNk_;
		MatrixOrIdentityType weRef(wftOptions_.twoSiteDmrg && ni>volumeOfNk_,we_);
		SizeType start = wsT_.getRowPtr(is);
		SizeType end = wsT_.getRowPtr(is+1);
		ComplexOrRealType sum=0;
		for (SizeType k2=weRef.getRowPtr(jen);k2<weRef.getRowPtr(jen+1);k2++) {
			int jpr = weRef.getColOrExit(k2);
			// jpr < 0 could be due to an m smaller than h, the Hilbert size of one site
			// this is checked against elsewhere
			assert(jpr >= 0);
			SizeType jp = dmrgWaveStruct_.lrs.right().permutationInverse(jpl + jpr*volumeOfNk_);
			ComplexOrRealType sum2 = 0;
			for (SizeType k = start;k < end;k++) {
				SizeType ip = wsT_.getCol(k);
				SizeType y = dmrgWaveStruct_.lrs.super().permutationInverse(ip + jp*nalpha);
				y -= offset;
				sum2 += wsT_.getValue(k)*src_.fastAccess(iOld_, y)*weRef.getValue(k2);
			}

			sum += sum2;
		}

		return sum;
	}

	ComplexOrRealType createAux1bFromInfinite(SizeType ip,
	                                          SizeType kp,
	                                          SizeType jp) const
	{
		SizeType offset = src_.offset(iOld_);
		const SparseMatrixType& ws = wsT_;
		const SparseMatrixType& weT = we_;
		SizeType ni= dmrgWaveStruct_.lrs.left().size();
		SizeType nip = dmrgWaveStruct_.lrs.left().permutationInverse().size()/volumeOfNk_;
		MatrixOrIdentityType wsRef2(wftOptions_.twoSiteDmrg && nip>volumeOfNk_, ws);
		SizeType start = weT.getRowPtr(jp);
		SizeType end = weT.getRowPtr(jp+1);
		ComplexOrRealType sum=0;
		for (SizeType k3=wsRef2.getRowPtr(ip);k3<wsRef2.getRowPtr(ip+1);k3++) {
			int ip2 = wsRef2.getColOrExit(k3);
			if (ip2<0) continue;
			SizeType alpha = dmrgWaveStruct_.lrs.left().permutationInverse(ip2+kp*nip);

			for (SizeType k = start; k < end; k++) {
				SizeType jp2 = weT.getCol(k);
				SizeType x = dmrgWaveStruct_.lrs.super().permutationInverse(alpha+jp2*ni);
				x -= offset;
				sum += weT.getValue(k)*src_.fastAccess(iOld_, x)*wsRef2.getValue(k3);
			}
		}

		return sum;
	}

	VectorWithOffsetType& dest_;
	SizeType i0_;
	const VectorWithOffsetType& src_;
	SizeType iOld_;
	const DmrgWaveStructType& dmrgWaveStruct_;
	const WftOptionsType& wftOptions_;
	const LeftRightSuperType& lrs_;
	const SparseMatrixType& wsT_;
	const SparseMatrixType& we_;
	SizeType volumeOfNk_;
	PackIndicesType pack1_;
	PackIndicesType pack2_;
	SizeType sysOrEnv_;

};
}
#endif // WFTSPARSETWOSITE_H
