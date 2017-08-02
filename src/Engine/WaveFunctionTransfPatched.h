/*
Copyright (c) 2009-2017, UT-Battelle, LLC
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

/*! \file WaveFunctionTransfPatched.h
 *
 *  This class implements the wave function transformation, see PRL 77, 3633 (1996)
 *  this is for when NOT using SU(2)
 *  and is optimized for patching
 */

#ifndef WFT_PATCHED_HEADER_H
#define WFT_PATCHED_HEADER_H

#include "PackIndices.h"
#include "ProgressIndicator.h"
#include "VectorWithOffsets.h" // so that PsimagLite::norm() becomes visible here
#include "VectorWithOffset.h" // so that PsimagLite::norm() becomes visible here
#include "WaveFunctionTransfBase.h"
#include "MatrixOrIdentity.h"
#include "MatrixVectorKron/GenIjPatch.h"

namespace Dmrg {

template<typename DmrgWaveStructType,typename VectorWithOffsetType>
class WaveFunctionTransfPatched : public
        WaveFunctionTransfBase<DmrgWaveStructType,VectorWithOffsetType> {

	typedef PsimagLite::PackIndices PackIndicesType;
	typedef WaveFunctionTransfBase<DmrgWaveStructType,VectorWithOffsetType> BaseType;
	typedef typename BaseType::VectorSizeType VectorSizeType;

public:

	typedef typename DmrgWaveStructType::BasisWithOperatorsType BasisWithOperatorsType;
	typedef typename BasisWithOperatorsType::SparseMatrixType SparseMatrixType;
	typedef typename BasisWithOperatorsType::BasisType BasisType;
	typedef typename SparseMatrixType::value_type SparseElementType;
	typedef typename PsimagLite::Vector<SparseElementType>::Type VectorType;
	typedef typename BasisWithOperatorsType::RealType RealType;
	typedef typename BasisType::FactorsType FactorsType;
	typedef typename DmrgWaveStructType::LeftRightSuperType LeftRightSuperType;
	typedef MatrixOrIdentity<SparseMatrixType> MatrixOrIdentityType;

	static const SizeType INFINITE = ProgramGlobals::INFINITE;
	static const SizeType EXPAND_SYSTEM = ProgramGlobals::EXPAND_SYSTEM;
	static const SizeType EXPAND_ENVIRON = ProgramGlobals::EXPAND_ENVIRON;

	WaveFunctionTransfPatched(const SizeType& stage,
	                        const bool& firstCall,
	                        const SizeType& counter,
	                        const DmrgWaveStructType& dmrgWaveStruct,
	                        bool twoSiteDmrg)
	    : stage_(stage),
	      firstCall_(firstCall),
	      counter_(counter),
	      dmrgWaveStruct_(dmrgWaveStruct),
	      twoSiteDmrg_(twoSiteDmrg),
	      progress_("WaveFunctionTransfPatched")
	{
		PsimagLite::OstringStream msg;
		msg<<"Constructing Local";
		progress_.printline(msg,std::cout);
	}

	virtual void transformVector(VectorWithOffsetType& psiDest,
	                             const VectorWithOffsetType& psiSrc,
	                             const LeftRightSuperType& lrs,
	                             const VectorSizeType& nk) const

	{
		if (stage_==EXPAND_ENVIRON) {
			if (firstCall_) {
				transformVector1FromInfinite(psiDest,psiSrc,lrs,nk);
			} else if (counter_==0) {
				transformVector1bounce(psiDest,psiSrc,lrs,nk);
			} else {
				transformVector1(psiDest,psiSrc,lrs,nk);
			}
		}

		if (stage_==EXPAND_SYSTEM) {
			if (firstCall_)
				transformVector2FromInfinite(psiDest,psiSrc,lrs,nk);
			else if (counter_==0)
				transformVector2bounce(psiDest,psiSrc,lrs,nk);
			else transformVector2(psiDest,psiSrc,lrs,nk);
		}
	}

private:

	template<typename SomeVectorType>
	void transformVector1(SomeVectorType& psiDest,
	                      const SomeVectorType& psiSrc,
	                      const LeftRightSuperType& lrs,
	                      const VectorSizeType& nk) const
	{
		if (twoSiteDmrg_)
			return transformVector1FromInfinite(psiDest,psiSrc,lrs,nk);

		typename ProgramGlobals::DirectionEnum dir1 = ProgramGlobals::EXPAND_ENVIRON;
		for (SizeType ii=0;ii<psiDest.sectors();ii++) {
			SizeType i0 = psiDest.sector(ii);
			transformVectorParallel(psiDest,psiSrc,lrs,i0,nk,dir1);
		}
	}

	template<typename SomeVectorType>
	void transformVectorParallel(SomeVectorType& psiDest,
	                             const SomeVectorType& psiSrc,
	                             const LeftRightSuperType& lrs,
	                             SizeType i0,
	                             const VectorSizeType& nk,
	                             typename ProgramGlobals::DirectionEnum dir) const
	{
		err("WaveFunctionTransfPatched: not implemented yet\n");
//		SizeType m = psiSrc.sector(i0);
//		SizeType target = dmrgWaveStruct_.lrs.super.qn(m);
//		GenIjPatch<LeftRightSuperType> genIjPatch(dmrgWaveStruct_.lrs, target);
//		VectorSizeType lv = genIjPatch(GenIjPatch<LeftRightSuperType>::LEFT);
//		VectorSizeType rv = genIjPatch(GenIjPatch<LeftRightSuperType>::RIGHT);
//		SizeType total = lv.size();
//		assert(total == rv.size());
		// if (dir_ == DmrgWaveStructType::DIR_2)
		// reshape psiSrc --> psiSrc(ip', beta')
		// Ws^T_(ip, ip') psiSrc(ip', beta') We^T(beta', beta)
		// (ip, beta) --> (ip, (kp, jp)) --> ((ip, kp), jp)  --> psiDest(alpha, jp)
		//
		// else
		// reshape psiSrc --> psiSrc(alpha', jp')
		// Ws(alpha, alpha') psiSrc(alpha', jp') We(jp', jp)
		// (alpha, jp) --> ((ip, kp), jp) --> (ip, (kp, jp)) --> psiDest(ip, beta)
	}

	template<typename SomeVectorType>
	void transformVector1FromInfinite(SomeVectorType& psiDest,
	                                  const SomeVectorType& psiSrc,
	                                  const LeftRightSuperType& lrs,
	                                  const VectorSizeType& nk) const
	{
		for (SizeType ii=0;ii<psiDest.sectors();ii++) {
			SizeType i0 = psiDest.sector(ii);
			transformVector1FromInfinite(psiDest,psiSrc,lrs,i0,nk);
		}
	}

	template<typename SomeVectorType>
	void transformVector1FromInfinite(SomeVectorType& psiDest,
	                                  const SomeVectorType& psiSrc,
	                                  const LeftRightSuperType& lrs,
	                                  SizeType i0,
	                                  const VectorSizeType& nk) const
	{
		SizeType volumeOfNk = DmrgWaveStructType::volumeOf(nk);
		SizeType nip = lrs.super().permutationInverse().size()/
		        lrs.right().permutationInverse().size();

		assert(lrs.left().permutationInverse().size()==volumeOfNk ||
		       lrs.left().permutationInverse().size()==dmrgWaveStruct_.ws.row());
		assert(lrs.right().permutationInverse().size()/volumeOfNk==dmrgWaveStruct_.we.col());

		SizeType start = psiDest.offset(i0);
		SizeType total = psiDest.effectiveSize(i0);

		SparseMatrixType ws(dmrgWaveStruct_.ws);
		SparseMatrixType we(dmrgWaveStruct_.we);
		SparseMatrixType weT;
		transposeConjugate(weT,we);

		PackIndicesType pack1(nip);
		PackIndicesType pack2(volumeOfNk);
		for (SizeType x=0;x<total;x++) {
			SizeType ip,beta,kp,jp;
			pack1.unpack(ip,beta,(SizeType)lrs.super().permutation(x+start));
			pack2.unpack(kp,jp,(SizeType)lrs.right().permutation(beta));
			psiDest.fastAccess(i0,x)=createAux1bFromInfinite(psiSrc,ip,kp,jp,ws,weT,nk);
		}
	}

	template<typename SomeVectorType>
	SparseElementType createAux1bFromInfinite(const SomeVectorType& psiSrc,
	                                          SizeType ip,
	                                          SizeType kp,
	                                          SizeType jp,
	                                          const SparseMatrixType& ws,
	                                          const SparseMatrixType& weT,
	                                          const VectorSizeType& nk) const
	{
		SizeType volumeOfNk = DmrgWaveStructType::volumeOf(nk);
		SizeType ni=dmrgWaveStruct_.lrs.left().size();
		SizeType nip = dmrgWaveStruct_.lrs.left().permutationInverse().size()/volumeOfNk;
		MatrixOrIdentityType wsRef2(twoSiteDmrg_ && nip>volumeOfNk,ws);
		SizeType start = weT.getRowPtr(jp);
		SizeType end = weT.getRowPtr(jp+1);
		SparseElementType sum=0;
		for (SizeType k3=wsRef2.getRowPtr(ip);k3<wsRef2.getRowPtr(ip+1);k3++) {
			int ip2 = wsRef2.getColOrExit(k3);
			if (ip2<0) continue;
			SizeType alpha = dmrgWaveStruct_.lrs.left().permutationInverse(ip2+kp*nip);

			for (SizeType k = start; k < end; k++) {
				SizeType jp2 = weT.getCol(k);
				SizeType x = dmrgWaveStruct_.lrs.super().permutationInverse(alpha+jp2*ni);
				sum += weT.getValue(k)*psiSrc.slowAccess(x)*wsRef2.getValue(k3);

			}
		}
		return sum;
	}

	template<typename SomeVectorType>
	void transformVector2(SomeVectorType& psiDest,
	                      const SomeVectorType& psiSrc,
	                      const LeftRightSuperType& lrs,
	                      const VectorSizeType& nk) const
	{
		if (twoSiteDmrg_)
			return transformVector2FromInfinite(psiDest,psiSrc,lrs,nk);

		typename ProgramGlobals::DirectionEnum dir2 = ProgramGlobals::EXPAND_SYSTEM;
		for (SizeType ii=0;ii<psiDest.sectors();ii++) {
			SizeType i0 = psiDest.sector(ii);
			transformVectorParallel(psiDest,psiSrc,lrs,i0,nk,dir2);
		}
	}

	template<typename SomeVectorType>
	void transformVector2FromInfinite(SomeVectorType& psiDest,
	                                  const SomeVectorType& psiSrc,
	                                  const LeftRightSuperType& lrs,
	                                  const VectorSizeType& nk) const
	{
		PsimagLite::OstringStream msg;
		msg<<" Destination sectors "<<psiDest.sectors();
		msg<<" Source sectors "<<psiSrc.sectors();
		progress_.printline(msg,std::cout);
		assert(dmrgWaveStruct_.lrs.super().permutationInverse().size()==psiSrc.size());

		VectorType psiV;
		for (SizeType srcI = 0; srcI < psiSrc.sectors(); ++srcI) {
			SizeType srcII = psiSrc.sector(srcI);
			psiSrc.extract(psiV,srcII);
			SizeType offset = psiSrc.offset(srcII);
			for (SizeType ii=0;ii<psiDest.sectors();ii++) {
				SizeType i0 = psiDest.sector(ii);
				SizeType start = psiDest.offset(i0);
				SizeType final = psiDest.effectiveSize(i0)+start;
				VectorType dest(final-start,0.0);
				if (srcI > 0) psiDest.extract(dest,i0);
				transformVector2FromInfinite(dest,start,psiV,offset,lrs,nk);
				psiDest.setDataInSector(dest,i0);
			}
		}
	}

	void transformVector2FromInfinite(VectorType& dest,
	                                  SizeType destOffset,
	                                  const VectorType& psiV,
	                                  SizeType offset,
	                                  const LeftRightSuperType& lrs,
	                                  const VectorSizeType& nk) const
	{
		SizeType volumeOfNk = DmrgWaveStructType::volumeOf(nk);
		SizeType nip = lrs.left().permutationInverse().size()/volumeOfNk;
		SizeType nalpha = lrs.left().permutationInverse().size();

		assert(nip==dmrgWaveStruct_.ws.col());

		const SparseMatrixType& we = dmrgWaveStruct_.we;
		const SparseMatrixType& ws = dmrgWaveStruct_.ws;
		SparseMatrixType wsT;
		transposeConjugate(wsT,ws);

		PackIndicesType pack1(nalpha);
		PackIndicesType pack2(nip);
		for (SizeType x=0;x<dest.size();x++) {
			SizeType isn,jen;
			pack1.unpack(isn,jen,(SizeType)lrs.super().permutation(x+destOffset));
			SizeType is,jpl;
			pack2.unpack(is,jpl,(SizeType)lrs.left().permutation(isn));

			dest[x] += createAux2bFromInfinite(psiV,offset,is,jpl,jen,wsT,we,nk);

		}
	}

	SparseElementType createAux2bFromInfinite(const VectorType& psiV,
	                                          SizeType offset,
	                                          SizeType is,
	                                          SizeType jpl,
	                                          SizeType jen,
	                                          const SparseMatrixType& wsT,
	                                          const SparseMatrixType& we,
	                                          const VectorSizeType& nk) const
	{
		SizeType nalpha=dmrgWaveStruct_.lrs.left().permutationInverse().size();
		SparseElementType sum=0;
		SizeType volumeOfNk = DmrgWaveStructType::volumeOf(nk);
		SizeType ni = dmrgWaveStruct_.lrs.right().size()/volumeOfNk;

		MatrixOrIdentityType weRef(twoSiteDmrg_ && ni>volumeOfNk,we);
		SizeType start = wsT.getRowPtr(is);
		SizeType end = wsT.getRowPtr(is+1);
		for (SizeType k2=weRef.getRowPtr(jen);k2<weRef.getRowPtr(jen+1);k2++) {
			int jpr = weRef.getColOrExit(k2);
			// jpr < 0 could be due to an m smaller than h, the Hilbert size of one site
			// this is checked against elsewhere
			assert(jpr >= 0);
			SizeType jp = dmrgWaveStruct_.lrs.right().permutationInverse(jpl + jpr*volumeOfNk);
			SparseElementType sum2 = 0;
			for (SizeType k = start;k < end;k++) {
				SizeType ip = wsT.getCol(k);
				SizeType y = dmrgWaveStruct_.lrs.super().permutationInverse(ip + jp*nalpha);
				if (y < offset || y - offset >= psiV.size()) continue;
				sum2 += wsT.getValue(k)*psiV[y-offset]*weRef.getValue(k2);
			}

			sum += sum2;
		}

		return sum;
	}

	template<typename SomeVectorType>
	void transformVector1bounce(SomeVectorType& psiDest,
	                            const SomeVectorType& psiSrc,
	                            const LeftRightSuperType& lrs,
	                            const VectorSizeType& nk) const
	{
		for (SizeType ii=0;ii<psiDest.sectors();ii++) {
			SizeType i0 = psiDest.sector(ii);
			transformVector1bounce(psiDest,psiSrc,lrs,i0,nk);
		}
	}

	template<typename SomeVectorType>
	void transformVector1bounce(SomeVectorType& psiDest,
	                            const SomeVectorType& psiSrc,
	                            const LeftRightSuperType& lrs,
	                            SizeType i0,
	                            const VectorSizeType& nk) const
	{
		SizeType volumeOfNk = DmrgWaveStructType::volumeOf(nk);
		SizeType nip = lrs.super().permutationInverse().size()/
		        lrs.right().permutationInverse().size();
		PsimagLite::OstringStream msg;
		msg<<" We're bouncing on the right, so buckle up!";
		progress_.printline(msg,std::cout);

		assert(dmrgWaveStruct_.lrs.super().permutationInverse().size()==psiSrc.size());

		SizeType start = psiDest.offset(i0);
		SizeType total = psiDest.effectiveSize(i0);

		SizeType nalpha=dmrgWaveStruct_.lrs.left().permutationInverse().size();
		PackIndicesType pack1(nip);
		PackIndicesType pack2(volumeOfNk);
		MatrixOrIdentityType wsRef(twoSiteDmrg_,dmrgWaveStruct_.ws);
		SizeType nip2 = (twoSiteDmrg_) ? dmrgWaveStruct_.ws.col() : nip;

		for (SizeType x=0;x<total;x++) {
			psiDest.fastAccess(i0,x) = 0.0;
			SizeType ip,beta,kp,jp;
			pack1.unpack(ip,beta,(SizeType)lrs.super().permutation(x+start));
			for (SizeType k=wsRef.getRowPtr(ip);k<wsRef.getRowPtr(ip+1);k++) {
				int ip2 = wsRef.getColOrExit(k);
				if (ip2 < 0) continue;
				pack2.unpack(kp,jp,(SizeType)lrs.right().permutation(beta));
				SizeType ipkp = dmrgWaveStruct_.lrs.left().permutationInverse(ip2 + kp*nip2);
				SizeType y = dmrgWaveStruct_.lrs.super().permutationInverse(ipkp + jp*nalpha);
				psiDest.fastAccess(i0,x) += psiSrc.slowAccess(y)*wsRef.getValue(k);
			}
		}
	}

	template<typename SomeVectorType>
	void transformVector2bounce(SomeVectorType& psiDest,
	                            const SomeVectorType& psiSrc,
	                            const LeftRightSuperType& lrs,
	                            const VectorSizeType& nk) const
	{
		for (SizeType ii=0;ii<psiDest.sectors();ii++) {
			SizeType i0 = psiDest.sector(ii);
			transformVector2bounce(psiDest,psiSrc,lrs,i0,nk);
		}
	}

	template<typename SomeVectorType>
	void transformVector2bounce(SomeVectorType& psiDest,
	                            const SomeVectorType& psiSrc,
	                            const LeftRightSuperType& lrs,
	                            SizeType i0,
	                            const VectorSizeType& nk) const
	{
		SizeType volumeOfNk = DmrgWaveStructType::volumeOf(nk);
		SizeType nip = lrs.left().permutationInverse().size()/volumeOfNk;
		SizeType nalpha = lrs.left().permutationInverse().size();

		PsimagLite::OstringStream msg;
		msg<<" We're bouncing on the left, so buckle up!";
		progress_.printline(msg,std::cout);

		assert(dmrgWaveStruct_.lrs.super().permutationInverse().size()==psiSrc.size());

		SizeType start = psiDest.offset(i0);
		SizeType total = psiDest.effectiveSize(i0);
		PackIndicesType pack1(nalpha);
		PackIndicesType pack2(nip);
		MatrixOrIdentityType weRef(twoSiteDmrg_,dmrgWaveStruct_.we);

		for (SizeType x=0;x<total;x++) {
			psiDest.fastAccess(i0,x) = 0.0;

			SizeType ip,alpha,kp,jp;
			pack1.unpack(alpha,jp,(SizeType)lrs.super().permutation(x+start));
			pack2.unpack(ip,kp,(SizeType)lrs.left().permutation(alpha));

			for (SizeType k=weRef.getRowPtr(jp);k<weRef.getRowPtr(jp+1);k++) {
				int jp2 = weRef.getColOrExit(k);
				if (jp2 < 0) continue;
				SizeType kpjp = dmrgWaveStruct_.lrs.right().
				        permutationInverse(kp + jp2*volumeOfNk);

				SizeType y = dmrgWaveStruct_.lrs.super().
				        permutationInverse(ip + kpjp*nip);
				psiDest.fastAccess(i0,x) += psiSrc.slowAccess(y) * weRef.getValue(k);
			}
		}
	}

	const SizeType& stage_;
	const bool& firstCall_;
	const SizeType& counter_;
	const DmrgWaveStructType& dmrgWaveStruct_;
	bool twoSiteDmrg_;
	PsimagLite::ProgressIndicator progress_;
}; // class WaveFunctionTransfPatched
} // namespace Dmrg

/*@}*/
#endif

