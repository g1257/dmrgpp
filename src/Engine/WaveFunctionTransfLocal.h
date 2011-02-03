// BEGIN LICENSE BLOCK
/*
Copyright (c) 2009, UT-Battelle, LLC
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
// END LICENSE BLOCK
/** \ingroup DMRG */
/*@{*/

/*! \file WaveFunctionTransformation.h
 *
 *  This class implements the wave function transformation, see PRL 77, 3633 (1996)
 *
 */

#ifndef WFT_LOCAL_HEADER_H
#define WFT_LOCAL_HEADER_H
 
#include "Utils.h"
#include "ProgressIndicator.h"
#include "VectorWithOffsets.h" // so that std::norm() becomes visible here
#include "VectorWithOffset.h" // so that std::norm() becomes visible here
#include "WaveFunctionTransfBase.h"

namespace Dmrg {
	

	template<typename DmrgWaveStructType,typename VectorWithOffsetType>
	class WaveFunctionTransfLocal : public
		WaveFunctionTransfBase<DmrgWaveStructType,VectorWithOffsetType> {
	public:
		typedef typename DmrgWaveStructType::BasisWithOperatorsType BasisWithOperatorsType;
		typedef typename BasisWithOperatorsType::SparseMatrixType SparseMatrixType;
		typedef typename BasisWithOperatorsType::BasisType BasisType;
		typedef typename SparseMatrixType::value_type SparseElementType;
		typedef std::vector<SparseElementType> VectorType;
		typedef typename BasisWithOperatorsType::RealType RealType;
		typedef typename BasisType::FactorsType FactorsType;
		static const size_t INFINITE = ProgramGlobals::INFINITE;
		static const size_t EXPAND_SYSTEM = ProgramGlobals::EXPAND_SYSTEM;
		static const size_t EXPAND_ENVIRON = ProgramGlobals::EXPAND_ENVIRON;
		
		WaveFunctionTransfLocal(
				const size_t& hilbertSpaceOneSite,
				const size_t& stage,
				const bool& firstCall,
				const size_t& counter,
				const DmrgWaveStructType& dmrgWaveStruct)
		: hilbertSpaceOneSite_(hilbertSpaceOneSite),
		  stage_(stage),
		  firstCall_(firstCall),
		  counter_(counter),
		  dmrgWaveStruct_(dmrgWaveStruct),
		  progress_("WaveFunctionTransfLocal",0)
		{}
		
		virtual void transformVector(
						VectorWithOffsetType& psiDest,
						const VectorWithOffsetType& psiSrc,
						const BasisWithOperatorsType& pSprime,
						const BasisWithOperatorsType& pEprime,
						const BasisType& pSE) const
		{
			//std::cerr<<"counter="<<counter_<<"direction = "<<stage_<<"\n";
			if (stage_==EXPAND_ENVIRON) {
				if (firstCall_) throw std::runtime_error("WFT: This corner case is unimplmemented yet (sorry!)\n");
				else if (counter_==0) transformVector1bounce(psiDest,psiSrc,pSprime,pEprime,pSE);
				else transformVector1(psiDest,psiSrc,pSprime,pEprime,pSE);
			}
			if (stage_==EXPAND_SYSTEM) {
				if (firstCall_) transformVector2FromInfinite(psiDest,psiSrc,pSprime,pEprime,pSE);
 				else if (counter_==0) transformVector2bounce(psiDest,psiSrc,pSprime,pEprime,pSE);
				else transformVector2(psiDest,psiSrc,pSprime,pEprime,pSE);
			}
		}

	private:
		
		template<typename SomeVectorType>
		void transformVector1(SomeVectorType& psiDest,const SomeVectorType& psiSrc,const BasisWithOperatorsType& pSprime,
				      const BasisWithOperatorsType& pEprime,const BasisType& pSE) const
		{
			for (size_t ii=0;ii<psiDest.sectors();ii++) {
				size_t i0 = psiDest.sector(ii);
				transformVector1(psiDest,psiSrc,pSprime,pEprime,pSE,i0);
			}
		}
		
		template<typename SomeVectorType>
		void transformVector1(SomeVectorType& psiDest,const SomeVectorType& psiSrc,const BasisWithOperatorsType& pSprime,
				      const BasisWithOperatorsType& pEprime,const BasisType& pSE,size_t i0) const
		{
			size_t nk = hilbertSpaceOneSite_;
			size_t nip = pSE.permutationInverse().size()/pEprime.permutationInverse().size();
			size_t njp = pEprime.permutationInverse().size()/nk;
			//printDmrgWave();
			if (dmrgWaveStruct_.pSprime.permutationInverse().size()!=dmrgWaveStruct_.ws.n_row()) {
				throw std::runtime_error("transformVector1():"
						"SpermutationInverse.size()!=dmrgWaveStruct_.ws.n_row()\n");
			}
			if (njp!=dmrgWaveStruct_.we.n_col()) {
				std::cerr<<"nip="<<nip<<" njp="<<njp<<" nk="<<nk<<" dmrgWaveStruct_.we.n_col()="<<dmrgWaveStruct_.we.n_col()<<"\n";
				throw std::runtime_error("WaveFunctionTransformation::transformVector1():"
						"njp!=dmrgWaveStruct_.we.n_col()\n");
			}
			/*if (dmrgWaveStruct_.pSE.permutationInverse().size()!=psiSrc.size()) {
				printDmrgWave();
				std::cerr<<"SEpermutationInverse.size="<<dmrgWaveStruct_.pSE.permutationInverse().size();
				std::cerr<<" psiSrc.size="<<psiSrc.size()<<"\n";
				throw std::runtime_error("WaveFunctionTransformation::transformVector1():"
						" dmrgWaveStruct_.SEpermutationInverse.size()!=dmrgWaveStruct_.psi.size()\n");
			}*/
			
			size_t start = psiDest.offset(i0);
			size_t final = psiDest.effectiveSize(i0)+start;
			
			SparseMatrixType ws(dmrgWaveStruct_.ws);
			SparseMatrixType we(dmrgWaveStruct_.we);
			SparseMatrixType weT;
			transposeConjugate(weT,we);
			
			for (size_t x=start;x<final;x++) {
				size_t ip,beta,kp,jp;
				utils::getCoordinates(ip,beta,(size_t)pSE.permutation(x),nip);
				utils::getCoordinates(kp,jp,(size_t)pEprime.permutation(beta),nk);
				psiDest[x]=createAux1b(psiSrc,ip,kp,jp,ws,weT);
			}
		}

		template<typename SomeVectorType>
		SparseElementType createAux1b(const SomeVectorType& psiSrc,size_t ip,size_t kp,size_t jp,
					const SparseMatrixType& ws,const SparseMatrixType& weT) const
		{
			size_t nk = hilbertSpaceOneSite_;
			size_t ni=dmrgWaveStruct_.ws.n_col();

			//int m = dmrgWaveStruct_.m;
			//size_t final = dmrgWaveStruct_.pSE.partition(m+1);
			//size_t start = dmrgWaveStruct_.pSE.partition(m);
			
			size_t nip = dmrgWaveStruct_.pSprime.permutationInverse().size()/nk;
			size_t alpha = dmrgWaveStruct_.pSprime.permutationInverse(ip+kp*nip);
			
			size_t offseta = 0; 
			size_t totala =dmrgWaveStruct_.ws.n_col();
			if (totala == dmrgWaveStruct_.pSprime.size()) {
				size_t eqn =  dmrgWaveStruct_.pSprime.qn(alpha,BasisType::BEFORE_TRANSFORM);
				int ma = dmrgWaveStruct_.pSprime.partitionFromQn(eqn);
				if (ma>=0) {
					offseta = dmrgWaveStruct_.pSprime.partition(ma);
					totala = dmrgWaveStruct_.pSprime.partition(ma+1) - offseta;
				}
			}

			size_t offsetb = 0;
			size_t totalb = dmrgWaveStruct_.pEprime.permutationInverse().size();
			if  (dmrgWaveStruct_.pEprime.dmrgTransformed()) {
				size_t eqn = dmrgWaveStruct_.pEprime.qn(jp);
				int mb =  dmrgWaveStruct_.pEprime.partitionFromQn(eqn,BasisType::BEFORE_TRANSFORM);
				if (mb>=0) {
					offsetb = dmrgWaveStruct_.pEprime.partition(mb,BasisType::BEFORE_TRANSFORM);
					totalb = dmrgWaveStruct_.pEprime.partition(mb+1,BasisType::BEFORE_TRANSFORM) - offsetb;
				}
			}
			
			SparseElementType sum=0;

			for (int k = ws.getRowPtr(alpha);k<ws.getRowPtr(alpha+1);k++) {
				size_t i = ws.getCol(k);
				for (int k2=weT.getRowPtr(jp);k2<weT.getRowPtr(jp+1);k2++) {
					size_t j = weT.getCol(k2);
					size_t x = dmrgWaveStruct_.pSE.permutationInverse(i+j*ni);
					sum += ws.getValue(k)*weT.getValue(k2)*psiSrc[x];
					//counter++;
				}
			}

			return sum;
		}
		
		template<typename SomeVectorType>
		void transformVector2(SomeVectorType& psiDest,const SomeVectorType& psiSrc,const BasisWithOperatorsType& pSprime,
				      const BasisWithOperatorsType& pEprime,const BasisType& pSE) const
		{
			for (size_t ii=0;ii<psiDest.sectors();ii++) {
				size_t i0 = psiDest.sector(ii);
				transformVector2(psiDest,psiSrc,pSprime,pEprime,pSE,i0);
			}
		}
		

		template<typename SomeVectorType>
		void transformVector2(SomeVectorType& psiDest,const SomeVectorType& psiSrc,const BasisWithOperatorsType& pSprime,
				      const BasisWithOperatorsType& pEprime,const BasisType& pSE,size_t i0) const
		{
			size_t nk = hilbertSpaceOneSite_;
			size_t nip = pSprime.permutationInverse().size()/nk;
			size_t nalpha = pSprime.permutationInverse().size();
			//printDmrgWave();
			if (dmrgWaveStruct_.pEprime.permutationInverse().size()!=dmrgWaveStruct_.we.n_row()) {
				throw std::runtime_error("transformVector2():"
						"PpermutationInverse.size()!=dmrgWaveStruct_.we.n_row()\n");
			}
			if (nip!=dmrgWaveStruct_.ws.n_col()) {
				throw std::runtime_error("WaveFunctionTransformation::transformVector2():"
						"nip!=dmrgWaveStruct_.ws.n_row()\n");
			}
			/*if (dmrgWaveStruct_.pSE.permutationInverse().size()!=psiSrc.size()) {
				printDmrgWave();
				std::cerr<<"SEpermutationInverse.size="<<dmrgWaveStruct_.pSE.permutationInverse().size();
				std::cerr<<" psiSrc.size="<<psiSrc.size()<<"\n";
				throw std::runtime_error("WaveFunctionTransformation::transformVector2():"
						" dmrgWaveStruct_.SEpermutationInverse.size()!=dmrgWaveStruct_.psi.size()\n");
			}*/

			size_t start = psiDest.offset(i0);
			size_t final = psiDest.effectiveSize(i0)+start;
			
			SparseMatrixType we(dmrgWaveStruct_.we);
			SparseMatrixType ws(dmrgWaveStruct_.ws);
			SparseMatrixType wsT;
			transposeConjugate(wsT,ws);
			
			for (size_t x=start;x<final;x++) {
				size_t ip,alpha,kp,jp;
				utils::getCoordinates(alpha,jp,(size_t)pSE.permutation(x),nalpha);
				utils::getCoordinates(ip,kp,(size_t)pSprime.permutation(alpha),nip);
				psiDest[x]=createAux2b(psiSrc,ip,kp,jp,wsT,we);
			}
			
		}

		template<typename SomeVectorType>
		SparseElementType createAux2b(const SomeVectorType& psiSrc,size_t ip,size_t kp,size_t jp,
					const SparseMatrixType& wsT,const SparseMatrixType& we) const
		{
			size_t nk = hilbertSpaceOneSite_;
			size_t nalpha=dmrgWaveStruct_.pSprime.permutationInverse().size();
			
			size_t beta = dmrgWaveStruct_.pEprime.permutationInverse(kp+jp*nk);
			
			SparseElementType sum=0;
			
			for (int k=wsT.getRowPtr(ip);k<wsT.getRowPtr(ip+1);k++) {
				//SparseElementType sum2=0;
				size_t alpha = wsT.getCol(k);
				for (int k2=we.getRowPtr(beta);k2<we.getRowPtr(beta+1);k2++) {
					size_t j = we.getCol(k2);
					size_t x = dmrgWaveStruct_.pSE.permutationInverse(alpha+j*nalpha);
					sum += wsT.getValue(k)*we.getValue(k2)*psiSrc[x];
				}
				//sum += sum2;
			}
			return sum;
		}
		
		template<typename SomeVectorType>
		void transformVector2FromInfinite(SomeVectorType& psiDest,const SomeVectorType& psiSrc,const BasisWithOperatorsType& pSprime,
				      const BasisWithOperatorsType& pEprime,const BasisType& pSE) const
		{
			for (size_t ii=0;ii<psiDest.sectors();ii++) {
				size_t i0 = psiDest.sector(ii);
				transformVector2FromInfinite(psiDest,psiSrc,pSprime,pEprime,pSE,i0);
			}
		}
		

		template<typename SomeVectorType>
		void transformVector2FromInfinite(SomeVectorType& psiDest,const SomeVectorType& psiSrc,const BasisWithOperatorsType& pSprime,
				      const BasisWithOperatorsType& pEprime,const BasisType& pSE,size_t i0) const
		{
			size_t nk = hilbertSpaceOneSite_;
			size_t nip = pSprime.permutationInverse().size()/nk;
			size_t nalpha = pSprime.permutationInverse().size();
			
			std::ostringstream msg;
			msg<<" We're moving to the finite loop, bumpy ride ahead!";
			progress_.printline(msg,std::cout);
			
			/*if (dmrgWaveStruct_.pEprime.permutationInverse().size()!=dmrgWaveStruct_.we.n_row()) {
				printDmrgWave();
				throw std::runtime_error("transformVector2():"
						"PpermutationInverse.size()!=dmrgWaveStruct_.we.n_row()\n");
			}*/
			if (nip!=dmrgWaveStruct_.ws.n_col()) {
				throw std::runtime_error("WaveFunctionTransformation::transformVector2():"
						"nip!=dmrgWaveStruct_.ws.n_row()\n");
			}
			if (dmrgWaveStruct_.pSE.permutationInverse().size()!=psiSrc.size()) {
				std::cerr<<"SEpermutationInverse.size="<<dmrgWaveStruct_.pSE.permutationInverse().size();
				std::cerr<<" psiSrc.size="<<psiSrc.size()<<"\n";
				throw std::runtime_error("WaveFunctionTransformation::transformVector2():"
						" dmrgWaveStruct_.SEpermutationInverse.size()!=dmrgWaveStruct_.psi.size()\n");
			}

			size_t start = psiDest.offset(i0);
			size_t final = psiDest.effectiveSize(i0)+start;
			
			SparseMatrixType we(dmrgWaveStruct_.we);
			SparseMatrixType ws(dmrgWaveStruct_.ws);
			SparseMatrixType wsT;
			transposeConjugate(wsT,ws);
			
			for (size_t x=start;x<final;x++) {
				size_t isn,jen;
				utils::getCoordinates(isn,jen,(size_t)pSE.permutation(x),nalpha);
				size_t is,jpl;
				utils::getCoordinates(is,jpl,(size_t)pSprime.permutation(isn),nip);
				//size_t jk,je;
				//utils::getCoordinates(jk,je,(size_t)pEprime.permutation(jen),npk);
				psiDest[x]=createAux2bFromInfinite(psiSrc,is,jpl,jen,wsT,we);
			}
			
		}
		
		
		// FIXME: INCOMING jen needs to be 4 times as big!!
		template<typename SomeVectorType>
		SparseElementType createAux2bFromInfinite(const SomeVectorType& psiSrc,size_t is,size_t jpl,size_t jen,
					const SparseMatrixType& wsT,const SparseMatrixType& we) const
		{
			size_t nk = hilbertSpaceOneSite_;
			size_t nalpha=dmrgWaveStruct_.pSprime.permutationInverse().size();
			
			
			
			SparseElementType sum=0;
			
			for (int k=wsT.getRowPtr(is);k<wsT.getRowPtr(is+1);k++) {
				//SparseElementType sum2=0;
				size_t ip = wsT.getCol(k);
				SparseElementType sum2 = 0;
				//for (int k2=we.getRowPtr(jen);k2<we.getRowPtr(jen+1);k2++) {
					size_t jpr = jen; //we.getCol(k2);
					size_t jp = dmrgWaveStruct_.pEprime.permutationInverse(jpl + jpr*nk);
					size_t y = dmrgWaveStruct_.pSE.permutationInverse(ip + jp*nalpha);
					sum2 += wsT.getValue(k)*psiSrc[y]; //*we.getValue(k2);
					
				//}
				sum += sum2;
			}
			return sum;
		}
		
		template<typename SomeVectorType>
		void transformVector1bounce(SomeVectorType& psiDest,const SomeVectorType& psiSrc,const BasisWithOperatorsType& pSprime,
				      const BasisWithOperatorsType& pEprime,const BasisType& pSE) const
		{
			for (size_t ii=0;ii<psiDest.sectors();ii++) {
				size_t i0 = psiDest.sector(ii);
				transformVector1bounce(psiDest,psiSrc,pSprime,pEprime,pSE,i0);
			}
		}
		
		template<typename SomeVectorType>
		void transformVector1bounce(SomeVectorType& psiDest,const SomeVectorType& psiSrc,const BasisWithOperatorsType& pSprime,
				      const BasisWithOperatorsType& pEprime,const BasisType& pSE,size_t i0) const
		{
			size_t nk = hilbertSpaceOneSite_;
			size_t nip = pSE.permutationInverse().size()/pEprime.permutationInverse().size();
			//size_t njp = pEprime.permutationInverse().size()/nk;
			//printDmrgWave();
			std::ostringstream msg;
			msg<<" We're bouncing on the right, so buckle up!";
			progress_.printline(msg,std::cout);
			
			if (dmrgWaveStruct_.pSE.permutationInverse().size()!=psiSrc.size()) {
				std::cerr<<"SEpermutationInverse.size="<<dmrgWaveStruct_.pSE.permutationInverse().size();
				std::cerr<<" psiSrc.size="<<psiSrc.size()<<"\n";
				throw std::runtime_error("WaveFunctionTransformation::transformVector1():"
						" dmrgWaveStruct_.SEpermutationInverse.size()!=dmrgWaveStruct_.psi.size()\n");
			}
			
			size_t start = psiDest.offset(i0);
			size_t final = psiDest.effectiveSize(i0)+start;
			
			/* SparseMatrixType ws(dmrgWaveStruct_.ws);
			SparseMatrixType we(dmrgWaveStruct_.we);
			SparseMatrixType weT;
			transposeConjugate(weT,we);
			*/
			//std::cerr<<"ALTERNATIVE\n";
			size_t nalpha=dmrgWaveStruct_.pSprime.permutationInverse().size();
			for (size_t x=start;x<final;x++) {
				size_t ip,beta,kp,jp;
				utils::getCoordinates(ip,beta,(size_t)pSE.permutation(x),nip);
				utils::getCoordinates(kp,jp,(size_t)pEprime.permutation(beta),nk);
				size_t ipkp = dmrgWaveStruct_.pSprime.permutationInverse(ip + kp*nip);
				size_t y = dmrgWaveStruct_.pSE.permutationInverse(ipkp + jp*nalpha);
				psiDest[x]=psiSrc[y];
			}
		}
		
		template<typename SomeVectorType>
		void transformVector2bounce(SomeVectorType& psiDest,const SomeVectorType& psiSrc,const BasisWithOperatorsType& pSprime,
				      const BasisWithOperatorsType& pEprime,const BasisType& pSE) const
		{
			for (size_t ii=0;ii<psiDest.sectors();ii++) {
				size_t i0 = psiDest.sector(ii);
				transformVector2bounce(psiDest,psiSrc,pSprime,pEprime,pSE,i0);
			}
		}
		

		template<typename SomeVectorType>
		void transformVector2bounce(SomeVectorType& psiDest,const SomeVectorType& psiSrc,const BasisWithOperatorsType& pSprime,
				      const BasisWithOperatorsType& pEprime,const BasisType& pSE,size_t i0) const
		{
			size_t nk = hilbertSpaceOneSite_;
			size_t nip = pSprime.permutationInverse().size()/nk;
			size_t nalpha = pSprime.permutationInverse().size();
			//printDmrgWave();
			
			std::ostringstream msg;
			msg<<" We're bouncing on the left, so buckle up!";
			progress_.printline(msg,std::cout);
			
			if (dmrgWaveStruct_.pSE.permutationInverse().size()!=psiSrc.size()) {
				std::cerr<<"SEpermutationInverse.size="<<dmrgWaveStruct_.pSE.permutationInverse().size();
				std::cerr<<" psiSrc.size="<<psiSrc.size()<<"\n";
				throw std::runtime_error("WaveFunctionTransformation::transformVector2():"
						" dmrgWaveStruct_.SEpermutationInverse.size()!=dmrgWaveStruct_.psi.size()\n");
			}

			size_t start = psiDest.offset(i0);
			size_t final = psiDest.effectiveSize(i0)+start;
			
			
			for (size_t x=start;x<final;x++) {
				size_t ip,alpha,kp,jp;
				utils::getCoordinates(alpha,jp,(size_t)pSE.permutation(x),nalpha);
				utils::getCoordinates(ip,kp,(size_t)pSprime.permutation(alpha),nip);
				size_t kpjp = dmrgWaveStruct_.pEprime.permutationInverse(kp + jp*nk);
				
				size_t y = dmrgWaveStruct_.pSE.permutationInverse(ip + kpjp*nip);
				psiDest[x]=psiSrc[y];
			}
			
		}

		const size_t& hilbertSpaceOneSite_;
		const size_t& stage_;
		const bool& firstCall_;
		const size_t& counter_;
		const DmrgWaveStructType& dmrgWaveStruct_;
		PsimagLite::ProgressIndicator progress_;
	}; // class WaveFunctionTransfLocal
} // namespace Dmrg

/*@}*/
#endif
