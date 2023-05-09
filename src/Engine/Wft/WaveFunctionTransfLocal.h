/*
Copyright (c) 2009-2013, 2017, UT-Battelle, LLC
All rights reserved

[DMRG++, Version 4.]
[by G.A., Oak Ridge National Laboratory]

UT Battelle Open Source Software License 11242008

OPEN SOURCE LICENSE

Subject to the conditions of this License, each
contributor to this software hereby grants, free of
charge, to any person obtaining a copy of this software
and associated documentation files (the "Software"), a
perpetual, worldwide, non-exclusive, no-charge,
royalty-free, irrevocable copyright license to use, copy,
modify, merge, publish, distribute, and/or sublicense
copies of the Software.

1. Redistributions of Software must retain the above
copyright and license notices, this list of conditions,
and the following disclaimer.  Changes or modifications
to, or derivative works of, the Software should be noted
with comments and the contributor and organization's
name.

2. Neither the names of UT-Battelle, LLC or the
Department of Energy nor the names of the Software
contributors may be used to endorse or promote products
derived from this software without specific prior written
permission of UT-Battelle.

3. The software and the end-user documentation included
with the redistribution, with or without modification,
must include the following acknowledgment:

"This product includes software produced by UT-Battelle,
LLC under Contract No. DE-AC05-00OR22725  with the
Department of Energy."

*********************************************************
DISCLAIMER

THE SOFTWARE IS SUPPLIED BY THE COPYRIGHT HOLDERS AND
CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED
WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A
PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE
COPYRIGHT OWNER, CONTRIBUTORS, UNITED STATES GOVERNMENT,
OR THE UNITED STATES DEPARTMENT OF ENERGY BE LIABLE FOR
ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE
OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH
DAMAGE.

NEITHER THE UNITED STATES GOVERNMENT, NOR THE UNITED
STATES DEPARTMENT OF ENERGY, NOR THE COPYRIGHT OWNER, NOR
ANY OF THEIR EMPLOYEES, REPRESENTS THAT THE USE OF ANY
INFORMATION, DATA, APPARATUS, PRODUCT, OR PROCESS
DISCLOSED WOULD NOT INFRINGE PRIVATELY OWNED RIGHTS.

*********************************************************

*/
/** \ingroup DMRG */
/*@{*/

/*! \file WaveFunctionTransfLocal.h
 *
 *  This class implements the wave function transformation, see PRL 77, 3633 (1996)
 *  this is for when NOT using SU(2)
 */

#ifndef WFT_LOCAL_HEADER_H
#define WFT_LOCAL_HEADER_H

#include "ProgressIndicator.h"
#include "VectorWithOffsets.h" // so that PsimagLite::norm() becomes visible here
#include "VectorWithOffset.h" // so that PsimagLite::norm() becomes visible here
#include "WaveFunctionTransfBase.h"
#include "MatrixOrIdentity.h"
#include "ParallelWftOne.h"
#include "Parallelizer.h"
#include "MatrixVectorKron/KronMatrix.h"
#include "WftAccelBlocks.h"
#include "WftAccelPatches.h"
#include "WftSparseTwoSite.h"
#include "WftAccelSvd.h"
#include "Profiling.h"

namespace Dmrg {

template<typename DmrgWaveStructType,
         typename VectorWithOffsetType,
         typename OptionsType_,
         typename OneSiteSpacesType_>
class WaveFunctionTransfLocal : public
        WaveFunctionTransfBase<DmrgWaveStructType,
        VectorWithOffsetType,
        OptionsType_,
        OneSiteSpacesType_> {


	typedef WaveFunctionTransfBase<DmrgWaveStructType,
	VectorWithOffsetType,
	OptionsType_,
	OneSiteSpacesType_>
	BaseType;
	typedef typename BaseType::VectorSizeType VectorSizeType;
	typedef typename BaseType::PackIndicesType PackIndicesType;

public:

	using OneSiteSpacesType = OneSiteSpacesType_;
	typedef typename BaseType::WftOptionsType WftOptionsType;
	typedef typename DmrgWaveStructType::BasisWithOperatorsType BasisWithOperatorsType;
	typedef typename BasisWithOperatorsType::SparseMatrixType SparseMatrixType;
	typedef typename BasisWithOperatorsType::BasisType BasisType;
	typedef typename BasisType::QnType QnType;
	typedef typename SparseMatrixType::value_type SparseElementType;
	typedef typename PsimagLite::Vector<SparseElementType>::Type VectorType;
	typedef typename BasisWithOperatorsType::RealType RealType;
	typedef typename DmrgWaveStructType::LeftRightSuperType LeftRightSuperType;
	typedef MatrixOrIdentity<SparseMatrixType> MatrixOrIdentityType;
	typedef ParallelWftOne<VectorWithOffsetType,
	DmrgWaveStructType,
	OneSiteSpacesType_> ParallelWftType;
	typedef PsimagLite::Matrix<SparseElementType> MatrixType;
	typedef WftAccelBlocks<BaseType> WftAccelBlocksType;
	typedef WftAccelPatches<BaseType> WftAccelPatchesType;
	typedef WftSparseTwoSite<BaseType, MatrixOrIdentityType> WftSparseTwoSiteType;
	typedef WftAccelSvd<BaseType> WftAccelSvdType;

	WaveFunctionTransfLocal(const DmrgWaveStructType& dmrgWaveStruct,
	                        const WftOptionsType& wftOptions)
	    : dmrgWaveStruct_(dmrgWaveStruct),
	      wftOptions_(wftOptions),
	      wftAccelBlocks_(dmrgWaveStruct, wftOptions),
	      wftAccelPatches_(dmrgWaveStruct, wftOptions),
	      wftAccelSvd_(dmrgWaveStruct, wftOptions)
	{}

	virtual void transformVector(VectorWithOffsetType& psiDest,
	                             const VectorWithOffsetType& psiSrc,
	                             const LeftRightSuperType& lrs,
	                             const OneSiteSpacesType& oneSiteSpaces) const

	{
		PsimagLite::Profiling profiling("WFT", std::cout);

		if (wftOptions_.dir == ProgramGlobals::DirectionEnum::EXPAND_ENVIRON) {
			if (wftOptions_.firstCall) {
				transformVector1FromInfinite(psiDest,psiSrc,lrs,oneSiteSpaces);
			} else if (wftOptions_.bounce) {
				transformVector1bounce(psiDest,psiSrc,lrs,oneSiteSpaces);
			} else {
				transformVector1(psiDest,psiSrc,lrs,oneSiteSpaces);
			}

			return;
		}

		if (wftOptions_.dir == ProgramGlobals::DirectionEnum::EXPAND_SYSTEM) {
			if (wftOptions_.firstCall)
				transformVector2FromInfinite(psiDest,psiSrc,lrs,oneSiteSpaces);
			else if (wftOptions_.bounce)
				transformVector2bounce(psiDest,psiSrc,lrs,oneSiteSpaces);
			else
				transformVector2(psiDest,psiSrc,lrs,oneSiteSpaces);

			return;
		}

		err("WFT Local: Stage is not EXPAND_ENVIRON or EXPAND_SYSTEM\n");
	}

private:

	void transformVector1(VectorWithOffsetType& psiDest,
	                      const VectorWithOffsetType& psiSrc,
	                      const LeftRightSuperType& lrs,
	                      const OneSiteSpacesType& oneSiteSpaces) const
	{
		if (wftOptions_.twoSiteDmrg)
			return transformVector1FromInfinite(psiDest,psiSrc,lrs,oneSiteSpaces);

		const ProgramGlobals::DirectionEnum dir1 = ProgramGlobals::DirectionEnum::EXPAND_ENVIRON;

		OneSiteSpacesType oneSiteSpaces2(oneSiteSpaces);
		oneSiteSpaces2.setDir(dir1);
		for (SizeType ii = 0; ii < psiDest.sectors(); ++ii)
			transformVectorParallel(psiDest, psiSrc, lrs, ii, oneSiteSpaces2);
	}

	void transformVectorParallel(VectorWithOffsetType& psiDest,
	                             const VectorWithOffsetType& psiSrc,
	                             const LeftRightSuperType& lrs,
	                             SizeType iNew,
	                             const OneSiteSpacesType& oneSiteSpaces) const
	{
		if (wftOptions_.accel == WftOptionsType::ACCEL_PATCHES) {
			SizeType iOld = findIold(psiSrc, psiDest.qn(iNew));
			return wftAccelPatches_(psiDest, iNew, psiSrc, iOld, lrs, oneSiteSpaces);
		}

		if (wftOptions_.accel == WftOptionsType::ACCEL_SVD) {
			SizeType iOld = findIold(psiSrc, psiDest.qn(iNew));
			return wftAccelSvd_(psiDest, iNew, psiSrc, iOld, lrs, oneSiteSpaces);
		}

		SizeType i0 = psiDest.sector(iNew);
		typedef PsimagLite::Parallelizer<ParallelWftType> ParallelizerType;

		ParallelizerType threadedWft(PsimagLite::Concurrency::codeSectionParams);
		ParallelWftType helperWft(psiDest,
		                          psiSrc,
		                          lrs,
		                          i0,
		                          oneSiteSpaces,
		                          dmrgWaveStruct_);

		threadedWft.loopCreate(helperWft);
	}

	void transformVector1FromInfinite(VectorWithOffsetType& psiDest,
	                                  const VectorWithOffsetType& psiSrc,
	                                  const LeftRightSuperType& lrs,
	                                  const OneSiteSpacesType& oneSiteSpaces) const
	{
		for (SizeType ii=0;ii<psiSrc.sectors();ii++) {
			SizeType iOld = psiSrc.sector(ii);
			const QnType& qn = psiSrc.qn(ii);
			SizeType iNew = psiDest.sector(findIold(psiDest, qn));
			tVector1FromInfinite(psiDest, iNew, psiSrc, iOld, lrs, oneSiteSpaces);
		}
	}

	void tVector1FromInfinite(VectorWithOffsetType& psiDest,
	                          SizeType i0,
	                          const VectorWithOffsetType& psiSrc,
	                          SizeType iOld,
	                          const LeftRightSuperType& lrs,
	                          const OneSiteSpacesType& oneSiteSpaces) const
	{
		if (wftOptions_.accel == WftOptionsType::ACCEL_BLOCKS &&
		        lrs.left().block().size() > 1)
			return wftAccelBlocks_.environFromInfinite(psiDest,
			                                           i0,
			                                           psiSrc,
			                                           iOld,
			                                           lrs,
			                                           oneSiteSpaces);

		typedef PsimagLite::Parallelizer<WftSparseTwoSiteType> ParallelizerType;

		SparseMatrixType ws;
		dmrgWaveStruct_.getTransform(ProgramGlobals::SysOrEnvEnum::SYSTEM).toSparse(ws);
		SparseMatrixType we;
		dmrgWaveStruct_.getTransform(ProgramGlobals::SysOrEnvEnum::ENVIRON).toSparse(we);
		SparseMatrixType weT;
		transposeConjugate(weT,we);

		ParallelizerType threadedWft(PsimagLite::Concurrency::codeSectionParams);

		WftSparseTwoSiteType helperWft(psiDest,
		                               i0,
		                               psiSrc,
		                               iOld,
		                               dmrgWaveStruct_,
		                               wftOptions_,
		                               lrs,
		                               oneSiteSpaces,
		                               ws,
		                               weT,
		                               ProgramGlobals::SysOrEnvEnum::ENVIRON);

		threadedWft.loopCreate(helperWft);
	}

	// --------------------------

	void transformVector2(VectorWithOffsetType& psiDest,
	                      const VectorWithOffsetType& psiSrc,
	                      const LeftRightSuperType& lrs,
	                      const OneSiteSpacesType& oneSiteSpaces) const
	{
		if (wftOptions_.twoSiteDmrg)
			return transformVector2FromInfinite(psiDest,psiSrc,lrs,oneSiteSpaces);

		const ProgramGlobals::DirectionEnum dir2 = ProgramGlobals::DirectionEnum::EXPAND_SYSTEM;
		OneSiteSpacesType oneSiteSpaces2(oneSiteSpaces);
		oneSiteSpaces2.setDir(dir2);

		for (SizeType ii=0;ii<psiDest.sectors();ii++)
			transformVectorParallel(psiDest, psiSrc, lrs, ii, oneSiteSpaces2);
	}

	void transformVector2FromInfinite(VectorWithOffsetType& psiDest,
	                                  const VectorWithOffsetType& psiSrc,
	                                  const LeftRightSuperType& lrs,
	                                  const OneSiteSpacesType& oneSiteSpaces) const
	{
		typedef PsimagLite::Parallelizer<WftSparseTwoSiteType> ParallelizerType;

		assert(dmrgWaveStruct_.lrs().super().permutationInverse().size() == psiSrc.size());

		bool inBlocks = (lrs.right().block().size() > 1 &&
		                 wftOptions_.accel == WftOptionsType::ACCEL_BLOCKS);
		SparseMatrixType we;
		dmrgWaveStruct_.getTransform(ProgramGlobals::SysOrEnvEnum::ENVIRON).toSparse(we);
		SparseMatrixType ws;
		dmrgWaveStruct_.getTransform(ProgramGlobals::SysOrEnvEnum::SYSTEM).toSparse(ws);
		SparseMatrixType wsT;
		transposeConjugate(wsT,ws);
		VectorType psiV;
		for (SizeType srcI = 0; srcI < psiSrc.sectors(); ++srcI) {
			SizeType srcII = psiSrc.sector(srcI);
			psiSrc.extract(psiV,srcII);
			for (SizeType ii=0;ii<psiDest.sectors();ii++) {
				SizeType i0 = psiDest.sector(ii);
				SizeType start = psiDest.offset(i0);
				SizeType final = psiDest.effectiveSize(i0)+start;
				VectorType dest(final-start,0.0);
				if (srcI > 0) psiDest.extract(dest,i0);
				if (inBlocks) {
					wftAccelBlocks_.systemFromInfinite(psiDest,
					                                   i0,
					                                   psiSrc,
					                                   srcII,
					                                   lrs,
					                                   oneSiteSpaces);
					continue;
				} else {
					ParallelizerType threadedWft(PsimagLite::Concurrency::codeSectionParams);

					WftSparseTwoSiteType helperWft(psiDest,
					                               i0,
					                               psiSrc,
					                               srcII,
					                               dmrgWaveStruct_,
					                               wftOptions_,
					                               lrs,
					                               oneSiteSpaces,
					                               wsT,
					                               we,
					                               ProgramGlobals::SysOrEnvEnum::SYSTEM);

					threadedWft.loopCreate(helperWft);
				}
			}
		}
	}

	// expand environ
	void transformVector1bounce(VectorWithOffsetType& psiDest,
	                            const VectorWithOffsetType& psiSrc,
	                            const LeftRightSuperType& lrs,
	                            const OneSiteSpacesType& oneSiteSpaces) const
	{
		SparseMatrixType ws;
		dmrgWaveStruct_.getTransform(ProgramGlobals::SysOrEnvEnum::SYSTEM).toSparse(ws);
		MatrixOrIdentityType wsRef(wftOptions_.twoSiteDmrg, ws);
		for (SizeType ii=0;ii<psiDest.sectors();ii++) {
			SizeType i0 = psiDest.sector(ii);
			tVector1bounce(psiDest,psiSrc,lrs,i0,oneSiteSpaces,wsRef);
		}
	}

	// expand environ
	void tVector1bounce(VectorWithOffsetType& psiDest,
	                    const VectorWithOffsetType& psiSrc,
	                    const LeftRightSuperType& lrs,
	                    SizeType i0,
	                    const OneSiteSpacesType& oneSiteSpaces,
	                    const MatrixOrIdentityType& wsRef) const
	{
		SizeType nip = lrs.super().permutationInverse().size()/
		        lrs.right().permutationInverse().size();

		assert(dmrgWaveStruct_.lrs().super().permutationInverse().size() == psiSrc.size());

		SizeType start = psiDest.offset(i0);
		SizeType total = psiDest.effectiveSize(i0);

		SizeType nalpha=dmrgWaveStruct_.lrs().left().permutationInverse().size();
		PackIndicesType pack1(nip);
		SizeType sizeOfHilbertForSiteAdded = oneSiteSpaces.hilbertMain();
		PackIndicesType pack2(sizeOfHilbertForSiteAdded);

		SizeType nip2 = (wftOptions_.twoSiteDmrg) ?
		            dmrgWaveStruct_.getTransform(ProgramGlobals::SysOrEnvEnum::SYSTEM).cols() :
		            nip;

		for (SizeType x=0;x<total;x++) {
			psiDest.fastAccess(i0,x) = 0.0;
			SizeType ip,beta,kp,jp;
			pack1.unpack(ip,beta,(SizeType)lrs.super().permutation(x+start));
			for (SizeType k=wsRef.getRowPtr(ip);k<wsRef.getRowPtr(ip+1);k++) {
				int ip2 = wsRef.getColOrExit(k);
				if (ip2 < 0) continue;
				pack2.unpack(kp,jp,(SizeType)lrs.right().permutation(beta));
				SizeType ipkp = dmrgWaveStruct_.lrs().left().permutationInverse(ip2 + kp*nip2);
				SizeType y = dmrgWaveStruct_.lrs().super().permutationInverse(ipkp + jp*nalpha);
				psiDest.fastAccess(i0,x) += psiSrc.slowAccess(y)*wsRef.getValue(k);
			}
		}
	}

	// expand system
	void transformVector2bounce(VectorWithOffsetType& psiDest,
	                            const VectorWithOffsetType& psiSrc,
	                            const LeftRightSuperType& lrs,
	                            const OneSiteSpacesType& oneSiteSpaces) const
	{
		SparseMatrixType we;
		dmrgWaveStruct_.getTransform(ProgramGlobals::SysOrEnvEnum::ENVIRON).toSparse(we);
		MatrixOrIdentityType weRef(wftOptions_.twoSiteDmrg, we);
		for (SizeType ii=0;ii<psiDest.sectors();ii++) {
			SizeType i0 = psiDest.sector(ii);
			tVector2bounce(psiDest,psiSrc,lrs,i0,oneSiteSpaces,weRef);
		}
	}

	// expand system
	void tVector2bounce(VectorWithOffsetType& psiDest,
	                    const VectorWithOffsetType& psiSrc,
	                    const LeftRightSuperType& lrs,
	                    SizeType i0,
	                    const OneSiteSpacesType& oneSiteSpaces,
	                    const MatrixOrIdentityType& weRef) const
	{
		SizeType sizeOfHilbertForSiteAdded = oneSiteSpaces.hilbertMain();
		SizeType nip = lrs.left().permutationInverse().size()/sizeOfHilbertForSiteAdded;
		SizeType nalpha = lrs.left().permutationInverse().size();

		assert(dmrgWaveStruct_.lrs().super().permutationInverse().size()==psiSrc.size());

		SizeType start = psiDest.offset(i0);
		SizeType total = psiDest.effectiveSize(i0);
		PackIndicesType pack1(nalpha);
		PackIndicesType pack2(nip);

		for (SizeType x=0;x<total;x++) {
			psiDest.fastAccess(i0,x) = 0.0;

			SizeType ip,alpha,kp,jp;
			pack1.unpack(alpha, jp, lrs.super().permutation(x+start));
			pack2.unpack(ip, kp, lrs.left().permutation(alpha));

			for (SizeType k=weRef.getRowPtr(jp);k<weRef.getRowPtr(jp+1);k++) {
				int jp2 = weRef.getColOrExit(k);
				if (jp2 < 0) continue;
				SizeType kpjp = dmrgWaveStruct_.lrs().right().
				        permutationInverse(kp + jp2*sizeOfHilbertForSiteAdded);

				SizeType y = dmrgWaveStruct_.lrs().super().
				        permutationInverse(ip + kpjp*nip);
				psiDest.fastAccess(i0,x) += psiSrc.slowAccess(y) * weRef.getValue(k);
			}
		}
	}

	static SizeType findIold(const VectorWithOffsetType& psiSrc,
	                         const QnType& qn)
	{
		SizeType sectors = psiSrc.sectors();
		for (SizeType i = 0; i < sectors; ++i)
			if (psiSrc.qn(i) == qn)
				return i;

		err("WaveFunctionTransfLocal::findIold(): Cannot find sector in old vector\n");
		throw PsimagLite::RuntimeError("UNREACHABLE\n");
	}

	const DmrgWaveStructType& dmrgWaveStruct_;
	const WftOptionsType& wftOptions_;
	WftAccelBlocksType wftAccelBlocks_;
	WftAccelPatchesType wftAccelPatches_;
	WftAccelSvdType wftAccelSvd_;
}; // class WaveFunctionTransfLocal
} // namespace Dmrg

/*@}*/
#endif

