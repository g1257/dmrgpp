/*
Copyright (c) 2009-2013, UT-Battelle, LLC
All rights reserved

[DMRG++, Version 2.0.0]
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

/*! \file WaveFunctionTransfSu2.h
 *
 *  This class implements the wave function transformation, see PRL 77, 3633 (1996)
 *
 */

#ifndef WFT_SU2_H
#define WFT_SU2_H

#include "PackIndices.h"
#include "ProgressIndicator.h"
#include "VectorWithOffsets.h" // so that std::norm() becomes visible here
#include "VectorWithOffset.h" // so that std::norm() becomes visible here
#include "WaveFunctionTransfBase.h"
#include "Random48.h"
#include "ParallelWftSu2.h"
#include "MatrixOrIdentity.h"

namespace Dmrg {

template<typename DmrgWaveStructType,typename VectorWithOffsetType>
class WaveFunctionTransfSu2  : public
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
	typedef ParallelWftSu2<VectorWithOffsetType,
	DmrgWaveStructType,
	LeftRightSuperType> ParallelWftType;
	typedef MatrixOrIdentity<SparseMatrixType> MatrixOrIdentityType;

	static const SizeType INFINITE = ProgramGlobals::INFINITE;
	static const SizeType EXPAND_SYSTEM = ProgramGlobals::EXPAND_SYSTEM;
	static const SizeType EXPAND_ENVIRON = ProgramGlobals::EXPAND_ENVIRON;

	WaveFunctionTransfSu2(
	        const SizeType& stage,
	        const bool& firstCall,
	        const SizeType& counter,
	        const DmrgWaveStructType& dmrgWaveStruct,
	        bool twoSiteDmrg)
	    : stage_(stage),
	      firstCall_(firstCall),
	      counter_(counter),
	      dmrgWaveStruct_(dmrgWaveStruct),
	      twoSiteDmrg_(twoSiteDmrg),
	      progress_("WaveFunctionTransfLocal")
	{
		PsimagLite::OstringStream msg;
		msg<<"Constructing SU(2)";
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

		RealType norm1 = std::norm(psiSrc);
		RealType norm2 = std::norm(psiDest);

		if (fabs(norm1-norm2)>1e-5) {
			PsimagLite::OstringStream msg;
			msg<<"WARNING: orig. norm= "<<norm1<<" resulting norm= "<<norm2;
			progress_.printline(msg,std::cout);
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

		typename ParallelWftType::DirectionEnum dir1 = ParallelWftType::DIR_1;
		for (SizeType ii=0;ii<psiDest.sectors();ii++) {
			SizeType i0 = psiDest.sector(ii);
			transformVectorParallel(psiDest,psiSrc,lrs,i0,nk,dir1);
		}
	}

	template<typename SomeVectorType>
	void transformVector2(SomeVectorType& psiDest,
	                      const SomeVectorType& psiSrc,
	                      const LeftRightSuperType& lrs,
	                      const VectorSizeType& nk) const
	{
		if (twoSiteDmrg_)
			return transformVector2FromInfinite(psiDest,psiSrc,lrs,nk);

		typename ParallelWftType::DirectionEnum dir2 = ParallelWftType::DIR_2;
		for (SizeType ii=0;ii<psiDest.sectors();ii++) {
			SizeType i0 = psiDest.sector(ii);
			transformVectorParallel(psiDest,psiSrc,lrs,i0,nk,dir2);
		}
	}

	template<typename SomeVectorType>
	void transformVectorParallel(SomeVectorType& psiDest,
	                             const SomeVectorType& psiSrc,
	                             const LeftRightSuperType& lrs,
	                             SizeType i0,
	                             const VectorSizeType& nk,
	                             typename ParallelWftType::DirectionEnum dir) const
	{
		SizeType total = psiDest.effectiveSize(i0);

		typedef PsimagLite::Parallelizer<ParallelWftType> ParallelizerType;
		ParallelizerType threadedWft(PsimagLite::Concurrency::npthreads,
		                             PsimagLite::MPI::COMM_WORLD);

		ParallelWftType helperWft(psiDest,
		                          psiSrc,
		                          lrs,
		                          i0,
		                          nk,
		                          dmrgWaveStruct_,
		                          dir);

		threadedWft.loopCreate(total, helperWft);
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
		SizeType volumeOfNk = ParallelWftType::volumeOf(nk);
		SizeType nip = lrs.super().permutationInverse().size()/
		        lrs.right().permutationInverse().size();

		assert(lrs.left().permutationInverse().size()==volumeOfNk ||
		       lrs.left().permutationInverse().size()==dmrgWaveStruct_.ws.row());
		assert(lrs.right().permutationInverse().size()/volumeOfNk==dmrgWaveStruct_.we.col());

		SizeType start = psiDest.offset(i0);
		SizeType total = psiDest.effectiveSize(i0);
		const FactorsType& factorsSE = lrs.super().getFactors();
		const FactorsType& factorsE = lrs.right().getFactors();

		SparseMatrixType factorsInverseSE, factorsInverseE;
		transposeConjugate(factorsInverseSE,factorsSE);
		transposeConjugate(factorsInverseE,factorsE);

		SparseMatrixType ws(dmrgWaveStruct_.ws);
		SparseMatrixType we(dmrgWaveStruct_.we);
		SparseMatrixType weT;
		transposeConjugate(weT,we);

		PackIndicesType pack1(nip);
		PackIndicesType pack2(volumeOfNk);
		for (SizeType x=0;x<total;x++) {
			psiDest.fastAccess(i0,x) = 0.0;
			SizeType ip,beta,kp,jp;
			SizeType xx = x + start;
			for (int kI = factorsInverseSE.getRowPtr(xx);
			     kI < factorsInverseSE.getRowPtr(xx+1);
			     kI++) {
				pack1.unpack(ip,beta,(SizeType)factorsInverseSE.getCol(kI));
				for (int k2I = factorsInverseE.getRowPtr(beta);
				     k2I < factorsInverseE.getRowPtr(beta+1);
				     k2I++) {
					pack2.unpack(kp,jp,(SizeType)factorsInverseE.getCol(k2I));
					psiDest.fastAccess(i0,x) += factorsInverseSE.getValue(kI)*
					        factorsInverseE.getValue(k2I)*
					        createAux1bFromInfinite(psiSrc,ip,kp,jp,ws,weT,nk);
				}
			}
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
		SizeType volumeOfNk = ParallelWftType::volumeOf(nk);
		SizeType ni=dmrgWaveStruct_.ws.col();
		SizeType nip = dmrgWaveStruct_.lrs.left().permutationInverse().size()/volumeOfNk;
		MatrixOrIdentityType wsRef2(twoSiteDmrg_ && nip>volumeOfNk,ws);
		const FactorsType& factorsS = dmrgWaveStruct_.lrs.left().getFactors();
		const FactorsType& factorsSE = dmrgWaveStruct_.lrs.super().getFactors();
		SparseElementType sum=0;

		for (SizeType k=wsRef2.getRowPtr(ip);k<wsRef2.getRowPtr(ip+1);k++) {
			int i = wsRef2.getColOrExit(k);
			if (i < 0) continue;
			SizeType ipkp=i+kp*nip;
			for (int k2I=factorsS.getRowPtr(ipkp);k2I<factorsS.getRowPtr(ipkp+1);k2I++) {
				SizeType alpha = factorsS.getCol(k2I);
				for (int k2=weT.getRowPtr(jp);k2<weT.getRowPtr(jp+1);k2++) {
					SizeType j = weT.getCol(k2);
					SizeType r = alpha+j*ni;
					for (int kI=factorsSE.getRowPtr(r);kI<factorsSE.getRowPtr(r+1);kI++) {
						SizeType x = factorsSE.getCol(kI);
						sum += wsRef2.getValue(k)*weT.getValue(k2)*psiSrc.slowAccess(x)*
						        factorsSE.getValue(kI)*factorsS.getValue(k2I);
					}
				}
			}
		}

		return sum;
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
	                                  SizeType start,
	                                  const VectorType& psiV,
	                                  SizeType offset,
	                                  const LeftRightSuperType& lrs,
	                                  const VectorSizeType& nk) const
	{
		SizeType volumeOfNk = ParallelWftType::volumeOf(nk);
		SizeType nip = lrs.left().permutationInverse().size()/volumeOfNk;
		SizeType nalpha = lrs.left().permutationInverse().size();

		assert(nip==dmrgWaveStruct_.ws.col());

		const FactorsType& factorsS = lrs.left().getFactors();
		const FactorsType& factorsSE = lrs.super().getFactors();
		SparseMatrixType factorsInverseSE, factorsInverseS;
		transposeConjugate(factorsInverseSE,factorsSE);
		transposeConjugate(factorsInverseS,factorsS);

		const SparseMatrixType& we = dmrgWaveStruct_.we;
		const SparseMatrixType& ws = dmrgWaveStruct_.ws;
		SparseMatrixType wsT;
		transposeConjugate(wsT,ws);

		PackIndicesType pack1(nalpha);
		PackIndicesType pack2(nip);
		for (SizeType x=0;x<dest.size();x++) {
			SizeType xx = x + offset;
			dest[x] = 0.0;
			SizeType isn,jen;
			for (int kI=factorsInverseSE.getRowPtr(xx);
			     kI<factorsInverseSE.getRowPtr(xx+1);
			     kI++) {
				pack1.unpack(isn,jen,(SizeType)factorsInverseSE.getCol(kI));
				SizeType is,jpl;
				for (int k2I=factorsInverseS.getRowPtr(isn);
				     k2I<factorsInverseS.getRowPtr(isn+1);
				     k2I++) {
					pack2.unpack(is,jpl,(SizeType)factorsInverseS.getCol(k2I));

					dest[x] += factorsInverseSE.getValue(kI)*
					        factorsInverseS.getValue(k2I)*
					        createAux2bFromInfinite(psiV,start,is,jpl,jen,wsT,we,nk);
				}
			}
		}
	}

	SparseElementType createAux2bFromInfinite(const VectorType& psiSrc,
	                                          SizeType start,
	                                          SizeType ip,
	                                          SizeType kp,
	                                          SizeType jp,
	                                          const SparseMatrixType& wsT,
	                                          const SparseMatrixType& we,
	                                          const VectorSizeType& nk) const
	{
		SizeType nalpha=dmrgWaveStruct_.lrs.left().permutationInverse().size();
		SparseElementType sum=0;
		SizeType volumeOfNk = ParallelWftType::volumeOf(nk);
		SizeType ni = dmrgWaveStruct_.lrs.right().size()/volumeOfNk;
		const FactorsType& factorsE = dmrgWaveStruct_.lrs.right().getFactors();
		const FactorsType& factorsSE = dmrgWaveStruct_.lrs.super().getFactors();
		MatrixOrIdentityType weRef(twoSiteDmrg_ && ni>volumeOfNk,we);
		SizeType end = start + psiSrc.size();
		SizeType kpjp = kp+jp*volumeOfNk;

		for (int k2I=factorsE.getRowPtr(kpjp);k2I<factorsE.getRowPtr(kpjp+1);k2I++) {
			SizeType beta = factorsE.getCol(k2I);
			for (int k=wsT.getRowPtr(ip);k<wsT.getRowPtr(ip+1);k++) {
				SizeType alpha = wsT.getCol(k);
				for (SizeType k2=weRef.getRowPtr(beta);k2<weRef.getRowPtr(beta+1);k2++) {
					int j = weRef.getColOrExit(k2);
					if (j < 0) continue;
					SizeType r = alpha + j*nalpha;
					for (int kI=factorsSE.getRowPtr(r);kI<factorsSE.getRowPtr(r+1);kI++) {
						SizeType x = factorsSE.getCol(kI);
						assert(x >= start && x < end);
						x -= start;
						sum += wsT.getValue(k)*weRef.getValue(k2)*psiSrc[x]*
						        factorsSE.getValue(kI)*factorsE.getValue(k2I);
					}
				}
			}
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
		SizeType volumeOfNk = ParallelWftType::volumeOf(nk);
		SizeType nip = lrs.super().permutationInverse().size()/
		        lrs.right().permutationInverse().size();

		assert(dmrgWaveStruct_.lrs.super().permutationInverse().size()==psiSrc.size());

		SizeType start = psiDest.offset(i0);
		SizeType total = psiDest.effectiveSize(i0);
		const FactorsType& factorsSE = lrs.super().getFactors();
		const FactorsType& factorsE = lrs.right().getFactors();

		SparseMatrixType factorsInverseSE, factorsInverseE;
		transposeConjugate(factorsInverseSE,factorsSE);
		transposeConjugate(factorsInverseE,factorsE);

		SparseMatrixType ws(dmrgWaveStruct_.ws);

		PackIndicesType pack1(nip);
		PackIndicesType pack2(volumeOfNk);
		for (SizeType x=0;x<total;x++) {
			psiDest.fastAccess(i0,x) = 0.0;
			SizeType ip,beta,kp,jp;
			SizeType xx = x + start;
			for (int kI = factorsInverseSE.getRowPtr(xx);
			     kI < factorsInverseSE.getRowPtr(xx+1);
			     kI++) {
				pack1.unpack(ip,beta,(SizeType)factorsInverseSE.getCol(kI));
				for (int k2I = factorsInverseE.getRowPtr(beta);
				     k2I < factorsInverseE.getRowPtr(beta+1);
				     k2I++) {
					pack2.unpack(kp,jp,(SizeType)factorsInverseE.getCol(k2I));
					psiDest.fastAccess(i0,x) += factorsInverseSE.getValue(kI)*
					        factorsInverseE.getValue(k2I)*
					        transformVector1bounceAux(psiSrc,ip,kp,jp,ws,nk);
				}
			}
		}
	}

	template<typename SomeVectorType>
	SparseElementType transformVector1bounceAux(const SomeVectorType& psiSrc,
	                                            SizeType ip,
	                                            SizeType kp,
	                                            SizeType jp,
	                                            const SparseMatrixType& ws,
	                                            const VectorSizeType& nk) const
	{
		SizeType volumeOfNk = ParallelWftType::volumeOf(nk);
		SizeType ni=dmrgWaveStruct_.ws.col();
		SizeType nip = dmrgWaveStruct_.lrs.left().permutationInverse().size()/volumeOfNk;
		MatrixOrIdentityType wsRef2(twoSiteDmrg_ && nip>volumeOfNk,ws);
		const FactorsType& factorsS = dmrgWaveStruct_.lrs.left().getFactors();
		const FactorsType& factorsSE = dmrgWaveStruct_.lrs.super().getFactors();
		SparseElementType sum=0;

		for (SizeType k=wsRef2.getRowPtr(ip);k<wsRef2.getRowPtr(ip+1);k++) {
			int i = wsRef2.getColOrExit(k);
			if (i < 0) continue;
			SizeType ipkp=i+kp*nip;
			for (int k2I=factorsS.getRowPtr(ipkp);k2I<factorsS.getRowPtr(ipkp+1);k2I++) {
				SizeType alpha = factorsS.getCol(k2I);
				SizeType j = jp;
				SizeType r = alpha+j*ni;
				for (int kI=factorsSE.getRowPtr(r);kI<factorsSE.getRowPtr(r+1);kI++) {
					SizeType x = factorsSE.getCol(kI);
					sum += wsRef2.getValue(k)*psiSrc.slowAccess(x)*
					        factorsSE.getValue(kI)*factorsS.getValue(k2I);
				}
			}
		}

		return sum;
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
		SizeType volumeOfNk = ParallelWftType::volumeOf(nk);
		SizeType nip = lrs.left().permutationInverse().size()/volumeOfNk;
		SizeType nalpha = lrs.left().permutationInverse().size();

		assert(dmrgWaveStruct_.lrs.super().permutationInverse().size()==psiSrc.size());

		const FactorsType& factorsS = lrs.left().getFactors();
		const FactorsType& factorsSE = lrs.super().getFactors();
		SparseMatrixType factorsInverseSE, factorsInverseS;
		transposeConjugate(factorsInverseSE,factorsSE);
		transposeConjugate(factorsInverseS,factorsS);
		SizeType start = psiDest.offset(i0);
		SizeType end = start + psiDest.effectiveSize(i0);
		const SparseMatrixType& we = dmrgWaveStruct_.we;

		PackIndicesType pack1(nalpha);
		PackIndicesType pack2(nip);
		for (SizeType x=start;x<end;x++) {
			psiDest.fastAccess(i0,x-start) = 0.0;
			SizeType isn,jen;
			for (int kI=factorsInverseSE.getRowPtr(x);
			     kI<factorsInverseSE.getRowPtr(x+1);
			     kI++) {
				pack1.unpack(isn,jen,(SizeType)factorsInverseSE.getCol(kI));
				SizeType is,jpl;
				for (int k2I=factorsInverseS.getRowPtr(isn);
				     k2I<factorsInverseS.getRowPtr(isn+1);
				     k2I++) {
					pack2.unpack(is,jpl,(SizeType)factorsInverseS.getCol(k2I));

					psiDest.fastAccess(i0,x-start) += factorsInverseSE.getValue(kI)*
					        factorsInverseS.getValue(k2I)*
					        transformVector2bounceAux(psiSrc,start,is,jpl,jen,we,nk);
				}
			}
		}
	}

	template<typename SomeVectorType>
	SparseElementType transformVector2bounceAux(const SomeVectorType& psiSrc,
	                                            SizeType start,
	                                            SizeType ip,
	                                            SizeType kp,
	                                            SizeType jp,
	                                            const SparseMatrixType& we,
	                                            const VectorSizeType& nk) const
	{
		SizeType nalpha=dmrgWaveStruct_.lrs.left().permutationInverse().size();
		SparseElementType sum=0;
		SizeType volumeOfNk = ParallelWftType::volumeOf(nk);
		SizeType ni = dmrgWaveStruct_.lrs.right().size()/volumeOfNk;
		const FactorsType& factorsE = dmrgWaveStruct_.lrs.right().getFactors();
		const FactorsType& factorsSE = dmrgWaveStruct_.lrs.super().getFactors();
		MatrixOrIdentityType weRef(twoSiteDmrg_ && ni>volumeOfNk,we);
		SizeType end = start + psiSrc.size();
		SizeType kpjp = kp+jp*volumeOfNk;

		for (int k2I=factorsE.getRowPtr(kpjp);k2I<factorsE.getRowPtr(kpjp+1);k2I++) {
			SizeType beta = factorsE.getCol(k2I);
			SizeType alpha = ip;
			for (SizeType k2=weRef.getRowPtr(beta);k2<weRef.getRowPtr(beta+1);k2++) {
				int j = weRef.getColOrExit(k2);
				if (j < 0) continue;
				SizeType r = alpha + j*nalpha;
				for (int kI=factorsSE.getRowPtr(r);kI<factorsSE.getRowPtr(r+1);kI++) {
					SizeType x = factorsSE.getCol(kI);
					assert(x >= start && x < end);
					sum += weRef.getValue(k2)*psiSrc.slowAccess(x)*
					        factorsSE.getValue(kI)*factorsE.getValue(k2I);
				}
			}
		}

		return sum;
	}

	const SizeType& stage_;
	const bool& firstCall_;
	const SizeType& counter_;
	const DmrgWaveStructType& dmrgWaveStruct_;
	bool twoSiteDmrg_;
	PsimagLite::ProgressIndicator progress_;

}; // class WaveFunctionTransfSu2
} // namespace Dmrg

/*@}*/
#endif

