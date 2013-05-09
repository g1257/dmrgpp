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

namespace Dmrg {
	

	template<typename DmrgWaveStructType,typename VectorWithOffsetType>
	class WaveFunctionTransfSu2  : public
		WaveFunctionTransfBase<DmrgWaveStructType,VectorWithOffsetType> {

		typedef PsimagLite::PackIndices PackIndicesType;

	public:
		typedef typename DmrgWaveStructType::BasisWithOperatorsType
			BasisWithOperatorsType;
		typedef typename BasisWithOperatorsType::SparseMatrixType
			SparseMatrixType;
		typedef typename BasisWithOperatorsType::BasisType BasisType;
		typedef typename SparseMatrixType::value_type SparseElementType;
		typedef typename PsimagLite::Vector<SparseElementType>::Type VectorType;
		typedef typename BasisWithOperatorsType::RealType RealType;
		typedef typename BasisType::FactorsType FactorsType;
		typedef typename DmrgWaveStructType::LeftRightSuperType
							LeftRightSuperType;

		static const size_t INFINITE = ProgramGlobals::INFINITE;
		static const size_t EXPAND_SYSTEM = ProgramGlobals::EXPAND_SYSTEM;
		static const size_t EXPAND_ENVIRON = ProgramGlobals::EXPAND_ENVIRON;

		WaveFunctionTransfSu2(
				const size_t& stage,
				const bool& firstCall,
				const size_t& counter,
				const DmrgWaveStructType& dmrgWaveStruct,
				bool twoSiteDmrg)
		: stage_(stage),
		  firstCall_(firstCall),
		  counter_(counter),
		  dmrgWaveStruct_(dmrgWaveStruct),
		  progress_("WaveFunctionTransfLocal",0)
		{
			if (twoSiteDmrg)
				throw PsimagLite::RuntimeError("SU(2) does not support two-site DMRG yet\n");
		}

		
		virtual void transformVector(VectorWithOffsetType& psiDest,
					     const VectorWithOffsetType& psiSrc,
					     const LeftRightSuperType& lrs,
					     const typename PsimagLite::Vector<size_t>::Type& nk) const
		{
			if (stage_==EXPAND_ENVIRON)
				transformVector1Su2(psiDest,psiSrc,lrs,nk);
			if (stage_==EXPAND_SYSTEM)
				transformVector2Su2(psiDest,psiSrc,lrs,nk);
		}

	private:
		
		template<typename SomeVectorType>
		void transformVector1Su2(SomeVectorType& psiDest,
					 const SomeVectorType& psiSrc,
					 const LeftRightSuperType& lrs,
					 const typename PsimagLite::Vector<size_t>::Type& nk) const
		{
			assert((size_t)dmrgWaveStruct_.lrs.super().getFactors().row()==psiSrc.size());

			for (size_t ii=0;ii<psiDest.sectors();ii++) {
				size_t i = psiDest.sector(ii);
				size_t start = psiDest.offset(i);
				size_t final = psiDest.effectiveSize(i)+start;
				transformVector1Su2(psiDest,psiSrc,lrs,start,final,nk);
			}
		}
		
		template<typename SomeVectorType>
		void transformVector1Su2(SomeVectorType& psiDest,
					 const SomeVectorType& psiSrc,
					 const LeftRightSuperType& lrs,
					 size_t start,
					 size_t final,
					 const typename PsimagLite::Vector<size_t>::Type& nk) const
		{
			size_t volumeOfNk = this->volumeOf(nk);
			const FactorsType& factorsSE = lrs.super().getFactors();
			const FactorsType& factorsSEOld = dmrgWaveStruct_.lrs.super().getFactors();
			const FactorsType& factorsE = lrs.right().getFactors();
			size_t nip = lrs.super().getFactors().row()/lrs.right().getFactors().row();
			
			FactorsType factorsInverseSE,factorsInverseE;
			transposeConjugate(factorsInverseSE,factorsSE);

			transposeConjugate(factorsInverseE,factorsE);
			
			SparseMatrixType ws(dmrgWaveStruct_.ws),we(dmrgWaveStruct_.we),weT;
			transposeConjugate(weT,we);
			
			PackIndicesType pack1(nip);
			PackIndicesType pack2(volumeOfNk);
			for (size_t x=start;x<final;x++) {
				psiDest[x] = 0;
				for (int kI = factorsInverseSE.getRowPtr(x);kI < factorsInverseSE.getRowPtr(x+1);kI++) {
					size_t ip,beta;
					pack1.unpack(ip,beta,(size_t)factorsInverseSE.getCol(kI));
					for (int k2I = factorsInverseE.getRowPtr(beta);k2I < factorsInverseE.getRowPtr(beta+1);k2I++) {
						size_t kp,jp;
						pack2.unpack(kp,jp,(size_t)factorsInverseE.getCol(k2I));
						psiDest[x] += createVectorAux1bSu2(psiSrc,ip,kp,jp,factorsSEOld,ws,weT,nk)*
								factorsInverseSE.getValue(kI)*factorsInverseE.getValue(k2I);
					}
				}
			}
		}
		
		template<typename SomeVectorType>
		SparseElementType createVectorAux1bSu2(const SomeVectorType& psiSrc,
						       size_t ip,
						       size_t kp,
						       size_t jp,
						       const FactorsType& factorsSE,
						       const SparseMatrixType& ws,
						       const SparseMatrixType& weT,
						       const typename PsimagLite::Vector<size_t>::Type& nk) const
		{
			size_t volumeOfNk = this->volumeOf(nk);
			size_t ni=dmrgWaveStruct_.ws.col();
			const FactorsType& factorsS = dmrgWaveStruct_.lrs.left().getFactors();
			SparseElementType sum=0;
			size_t nip = dmrgWaveStruct_.lrs.left().permutationInverse().size()/volumeOfNk;
			
			size_t ipkp=ip+kp*nip;
			for (int k2I=factorsS.getRowPtr(ipkp);k2I<factorsS.getRowPtr(ipkp+1);k2I++) {
				size_t alpha = factorsS.getCol(k2I);
				for (int k=ws.getRowPtr(alpha);k<ws.getRowPtr(alpha+1);k++) {
					size_t i = ws.getCol(k);
					for (int k2=weT.getRowPtr(jp);k2<weT.getRowPtr(jp+1);k2++) {
						size_t j = weT.getCol(k2);
						size_t r = i+j*ni;
						for (int kI=factorsSE.getRowPtr(r);kI<factorsSE.getRowPtr(r+1);kI++) {
							size_t x = factorsSE.getCol(kI);
							sum += ws.getValue(k)*weT.getValue(k2)*psiSrc[x]*
								factorsSE.getValue(kI)*factorsS.getValue(k2I);
						}
					}
				}
			}
			return sum;
		}
		
		template<typename SomeVectorType>
		void transformVector2Su2(SomeVectorType& psiDest,
					 const SomeVectorType& psiSrc,
					 const LeftRightSuperType& lrs,
					 const typename PsimagLite::Vector<size_t>::Type& nk) const
		{
			assert(dmrgWaveStruct_.lrs.super().permutationInverse().size()==psiSrc.size());

			for (size_t ii=0;ii<psiDest.sectors();ii++) {
				size_t i = psiDest.sector(ii);
				size_t start = psiDest.offset(i);
				size_t final = psiDest.effectiveSize(i)+start;
				transformVector2Su2(psiDest,psiSrc,lrs,start,final,nk);
			}
		}
		
		template<typename SomeVectorType>
		void transformVector2Su2(SomeVectorType& psiDest,
					 const SomeVectorType& psiSrc,
					 const LeftRightSuperType& lrs,
					 size_t start,
					 size_t final,
					 const typename PsimagLite::Vector<size_t>::Type& nk) const
		{
			size_t volumeOfNk = this->volumeOf(nk);
			size_t nip = lrs.left().getFactors().row()/volumeOfNk;
			size_t nalpha = lrs.left().getFactors().row();
			
			const FactorsType& factorsSE = lrs.super().getFactors();

			const FactorsType& factorsS = lrs.left().getFactors();
			FactorsType factorsInverseSE,factorsInverseS;
			transposeConjugate(factorsInverseSE,factorsSE);

			transposeConjugate(factorsInverseS,factorsS);
			SparseMatrixType ws(dmrgWaveStruct_.ws),we(dmrgWaveStruct_.we),wsT;
			transposeConjugate(wsT,ws);
			if (dmrgWaveStruct_.lrs.left().getFactors().row()==0) {
				wsT.makeDiagonal(volumeOfNk,1.0);
			}
			PackIndicesType pack1(nalpha);
			PackIndicesType pack2(nip);
			for (size_t x=start;x<final;x++) {
				psiDest[x] = 0;
				size_t xx = x;
				for (int kI=factorsInverseSE.getRowPtr(xx);kI<factorsInverseSE.getRowPtr(xx+1);kI++) {
					size_t alpha,jp;
					pack1.unpack(alpha,jp,(size_t)factorsInverseSE.getCol(kI));
					size_t alphax =  alpha;
					for (int k2I=factorsInverseS.getRowPtr(alphax);k2I<factorsInverseS.getRowPtr(alphax+1);k2I++) {
						size_t ip,kp;
						pack2.unpack(ip,kp,(size_t)factorsInverseS.getCol(k2I));
						psiDest[x] += fastAux2bSu2(psiSrc,ip,kp,jp,wsT,we,nk)*
								factorsInverseSE.getValue(kI)*factorsInverseS.getValue(k2I);
					}
				}
			}
			RealType x = std::norm(psiDest);
			if (fabs(x)>1e-6) return;
//			assert(false);
			PsimagLite::OstringStream msg1,msg2,msg3;
			msg1<<"WARNING: WFT forced to stop, this is a know issue with SU(2).";
			progress_.printline(msg1,std::cout);
			msg2<<"WARNING: If you are targeting anything other the ground state results will be wrong";
			progress_.printline(msg2,std::cout);
			// hack to make su(2) work with the wft
			// Note that this will ``drop'' vectors on the floor
			// and will cause not trivial targeting  (like tst, dynamic, correction vector) to fail
			PsimagLite::Random48<RealType> random(34343);
			for (size_t x=start;x<final;x++) {
				psiDest[x] = random();
			}
			msg3<<"WARNING: I have dropped your vector on the floor!";
			progress_.printline(msg3,std::cout);
		}
		
		template<typename SomeVectorType>
		SparseElementType fastAux2bSu2(const SomeVectorType& psiSrc,
					       size_t ip,
					       size_t kp,
					       size_t jp,
					       const SparseMatrixType& wsT,
					       const SparseMatrixType& we,
					       const typename PsimagLite::Vector<size_t>::Type& nk) const
		{
			size_t nalpha=wsT.row(); //dmrgWaveStruct_.lrs.left().getFactors().rank();
			assert(nalpha>0);
			SparseElementType sum=0;
			const FactorsType& factorsE = dmrgWaveStruct_.lrs.right().getFactors();
			const FactorsType& factorsSE = dmrgWaveStruct_.lrs.super().getFactors();
			size_t volumeOfNk = this->volumeOf(nk);
			size_t kpjp = kp+jp*volumeOfNk;
			assert(kpjp<dmrgWaveStruct_.lrs.right().permutationInverse().size());
			size_t kpjpx = dmrgWaveStruct_.lrs.right().permutationInverse(kpjp);
			
			for (int k2I=factorsE.getRowPtr(kpjpx);k2I<factorsE.getRowPtr(kpjpx+1);k2I++) {
				size_t beta = factorsE.getCol(k2I);
				for (int k=wsT.getRowPtr(ip);k<wsT.getRowPtr(ip+1);k++) {
					size_t alpha = wsT.getCol(k);
					for (int k2=we.getRowPtr(beta);k2<we.getRowPtr(beta+1);k2++) {
						size_t j = we.getCol(k2);
						size_t r = alpha + j*nalpha;
						for (int kI=factorsSE.getRowPtr(r);kI<factorsSE.getRowPtr(r+1);kI++) {
							size_t x = factorsSE.getCol(kI);
							sum += wsT.getValue(k)*we.getValue(k2)*psiSrc[x]*
								factorsSE.getValue(kI)*factorsE.getValue(k2I);
						}
					}
				}
			}
			return sum;
		}

		const size_t& stage_;
		const bool& firstCall_;
		const size_t& counter_;
		const DmrgWaveStructType& dmrgWaveStruct_;
		PsimagLite::ProgressIndicator progress_;

	}; // class WaveFunctionTransfSu2
} // namespace Dmrg

/*@}*/
#endif
