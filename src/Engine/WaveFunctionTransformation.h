// BEGIN LICENSE BLOCK
/*
Copyright © 2009 , UT-Battelle, LLC
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

#ifndef WAVEFUNCTIONT_HEADER_H
#define WAVEFUNCTIONT_HEADER_H
 
#include "Utils.h"
#include "ProgressIndicator.h"
#include "VectorWithOffsets.h" // so that std::norm() becomes visible here
#include "VectorWithOffset.h" // so that std::norm() becomes visible here

namespace Dmrg {
	
	template<typename BasisWithOperatorsType>
	struct DmrgWaveStructure {
		typedef typename BasisWithOperatorsType::OperatorType OperatorType;
		typedef typename OperatorType::SparseMatrixType SparseMatrixType;
		typedef typename SparseMatrixType::value_type SparseElementType;
		
		psimag::Matrix<SparseElementType> ws;
		psimag::Matrix<SparseElementType> we;
		typename BasisWithOperatorsType::BasisType pSE;
		BasisWithOperatorsType pSprime,pEprime;
		//int m;
		//std::vector<SparseElementType> psi;

		DmrgWaveStructure() : pSE("pSE"),pSprime("pSprime"),pEprime("pEprime") { }
	};

	template<typename BasisWithOperatorsType>
	std::ostream& operator<<(std::ostream& os,
			const DmrgWaveStructure<BasisWithOperatorsType>& dmrgWave)
	{
		os<<"ws.nrow="<<dmrgWave.ws.n_row()<<" ws.ncol="<<dmrgWave.ws.n_col()<<"\n";
		os<<"we.nrow="<<dmrgWave.we.n_row()<<" we.ncol="<<dmrgWave.we.n_col()<<"\n";
		os<<"SpermutationInverse.size="<<dmrgWave.pSprime.permutationInverse().size()<<"\n";
		os<<"EpermutationInverse.size="<<dmrgWave.pEprime.permutationInverse().size()<<"\n";
		os<<"SEpermutationInverse.size="<<dmrgWave.pSE.permutationInverse().size()<<"\n";
		//os<<"psi.size="<<dmrgWave.psi.size()<<"\n";
		
		return os;
	}

	template<typename BasisWithOperatorsType>
	class WaveFunctionTransformation {
		public:
		
		enum {INFINITE=0,SHRINK_SYSTEM=1,SHRINK_ENVIRON=2};
		enum {DO_NOT_RESET_COUNTER,RESET_COUNTER};

		typedef typename BasisWithOperatorsType::SparseMatrixType SparseMatrixType;
		typedef typename BasisWithOperatorsType::BasisType BasisType;
		typedef typename SparseMatrixType::value_type SparseElementType;
		typedef std::vector<SparseElementType> VectorType;
		typedef typename BasisWithOperatorsType::RealType RealType;
		typedef typename BasisType::FactorsType FactorsType;
		
		WaveFunctionTransformation() 
			: isEnabled_(true),useSu2Symmetry_(BasisType::useSu2Symmetry()),
				     stage_(INFINITE),counter_(0),doNextOne_(true),firstCall_(true),progress_("WaveFunctionTransformation",0)
		{}

		void init(size_t nk)
		{
			sizeOfOneSiteHilbertSpace_=nk;
			useSu2Symmetry_ = BasisType::useSu2Symmetry();
		}

		void setStage(int stage,int option=RESET_COUNTER)
		{
			stage_=stage;
			if (option==DO_NOT_RESET_COUNTER) return;
			counter_=0;
		}
		
		void triggerOn(const BasisWithOperatorsType& pSprime,
	  			const BasisWithOperatorsType& pEprime,
      				const BasisType& pSE)
		{
			bool allow=false;
			switch (stage_) {
				case INFINITE:
					allow=false;
					break;
				case SHRINK_ENVIRON:
					allow=true;
					
				case SHRINK_SYSTEM:
					allow=true;
			}
			// FIXME: Must check the below change when using SU(2)!!
			//if (m<0) allow = false; // isEnabled_=false;
			
			if (!isEnabled_ || !allow) return;
			try {
				beforeWft(pSprime,pEprime,pSE);
			} catch (std::exception& e) {
				doNextOne_=false;
			}
			std::ostringstream msg;
			msg<<"I'm ready to take passengers";
			progress_.printline(msg,std::cout);
		}
		
		// FIXME: change name to transformVector
		template<typename SomeVectorType,typename SomeVectorType2>
		void setInitialVector(	
					SomeVectorType& dest,
					const SomeVectorType2& src,
					const BasisWithOperatorsType& pSprime,
	  				const BasisWithOperatorsType& pEprime,
      					const BasisType& pSE) const
		{
			bool allow=false;
			switch (stage_) {
				case INFINITE:
					allow=false;
					break;
				case SHRINK_ENVIRON:
					allow=true;
					
				case SHRINK_SYSTEM:
					allow=true;
			}
			// FIXME: Must check the below change when using SU(2)!!
			//if (m<0) allow = false; // isEnabled_=false;
			
			if (isEnabled_ && allow) {
				//try {
				//	beforeWft(pSprime,pEprime,pSE,m);
				//} catch (std::exception& e) {
				//	afterWft(pSprime,pEprime,pSE,m);
				
				if (!doNextOne_) {
					createRandomVector(dest);
					return;
				}
				RealType eps = 1e-6;
				if (std::norm(src)<eps) throw std::runtime_error("src's norm is zero\n");
				createVector(dest,src,pSprime,pEprime,pSE);
				if (std::norm(dest)<eps) throw std::runtime_error("dest's norm is zero\n");
				//afterWft(pSprime,pEprime,pSE,m);	
			} else {
				createRandomVector(dest);
			}
		}
		
		void triggerOff(const BasisWithOperatorsType& pSprime,
	  			const BasisWithOperatorsType& pEprime,
      				const BasisType& pSE) //,int m)
		{
			bool allow=false;
			switch (stage_) {
				case INFINITE:
					allow=false;
					break;
				case SHRINK_ENVIRON:
					allow=true;
					
				case SHRINK_SYSTEM:
					allow=true;
			}
			// FIXME: Must check the below change when using SU(2)!!
			//if (m<0) allow = false; // isEnabled_=false;
			
			if (!isEnabled_ || !allow) return;
			afterWft(pSprime,pEprime,pSE); //,m);	
			doNextOne_=true;
			std::ostringstream msg;
			msg<<"No more passengers, please";
			progress_.printline(msg,std::cout);
		}
		
		template<typename SomeVectorType>
		void createRandomVector(SomeVectorType& y) const
		{
			for (size_t jj=0;jj<y.sectors();jj++) {
				size_t j = y.sector(jj);
				size_t offset = y.offset(j);
				size_t total = y.effectiveSize(j);
				size_t final = offset + total;
				createRandomVector(y,offset,final);
			}
			if (!isEnabled_) return; // don't make noise unless enabled
			std::ostringstream msg;
			msg<<"Yes, I'm awake, but there's nothing heavy to do now";
			progress_.printline(msg,std::cout);
		}
		
		template<typename SomeVectorType>
		void createRandomVector(SomeVectorType& y,size_t offset,size_t final) const
		{
			typename SomeVectorType::value_type tmp;
   			RealType atmp=0;
			for (size_t i=offset;i<final;i++) {
				utils::myRandomT(tmp);
				y[i]=tmp;
				atmp += utils::myProductT(y[i],y[i]);
			}
			atmp = 1.0 / sqrt (atmp);
			for (size_t i=offset;i<final;i++) y[i] *= atmp;
		}

		template<typename SomeMatrixType>
		void push(
			const SomeMatrixType& transform,
			int option,
   			//const SomeVectorType& psi,
      			const BasisWithOperatorsType& pBasis,
			const BasisWithOperatorsType& pBasisSummed,
   			const BasisType& pSE)
//			 size_t m)
		{
			if (!isEnabled_) return;
			
			switch (stage_) {
				case INFINITE:
					if (option==0) {
						wsStack_.push(transform);
						dmrgWave_.ws=transform;
					} else {
						weStack_.push(transform);
						dmrgWave_.we=transform;
						//std::cerr<<"CHANGED dmrgWave_.we to transform\n";
						//std::cerr<<"PUSHING "<<transform.n_row()<<"x"<<transform.n_col()<<"\n";
					}
					break;
				case SHRINK_SYSTEM:
					if (option==0) throw std::logic_error("SHRINK_SYSTEM but option==0\n");
					dmrgWave_.we=transform;
					dmrgWave_.ws=transform;
					//vectorConvert(dmrgWave_.psi,psi);
					weStack_.push(transform);
					//std::cerr<<"PUSHING (POPPING) We "<<weStack_.size()<<"\n";
					break;
				case SHRINK_ENVIRON:
					if (option==1) throw std::logic_error("SHRINK_ENVIRON but option==1\n");
					dmrgWave_.ws=transform;
					dmrgWave_.we=transform;
					//vectorConvert(dmrgWave_.psi,psi);
					wsStack_.push(transform);
					break;
			}

			dmrgWave_.pSE=pSE;
			if (option==0) { // transforming the system
				dmrgWave_.pEprime=pBasisSummed;
				dmrgWave_.pSprime=pBasis;
			} else {
				dmrgWave_.pSprime=pBasisSummed;
				dmrgWave_.pEprime=pBasis;
			}
			std::ostringstream msg;
			msg<<"OK, pushing option="<<option<<" and stage="<<stage_;
			progress_.printline(msg,std::cout);
		}

		void disable() { isEnabled_=false; }
		
		bool isEnabled() const { return isEnabled_; }
		
	private:			
		bool isEnabled_;
		bool useSu2Symmetry_;
		size_t sizeOfOneSiteHilbertSpace_;
		size_t stage_;
		size_t counter_;
		bool doNextOne_;
		bool firstCall_;
		ProgressIndicator progress_;
		DmrgWaveStructure<BasisWithOperatorsType> dmrgWave_;
		std::stack<psimag::Matrix<SparseElementType> > wsStack_,weStack_;
		
		void beforeWft(const BasisWithOperatorsType& pSprime,
				  const BasisWithOperatorsType& pEprime,const BasisType& pSE)
		{
			if (stage_==SHRINK_SYSTEM) {
				if (wsStack_.size()>=1) {
					dmrgWave_.ws=wsStack_.top();
					wsStack_.pop();
				} else {
					std::cerr<<"PUSHING STACK ERROR S\n";
					throw std::runtime_error("System Stack is empty\n");
				}
			}
			
			if (stage_==SHRINK_ENVIRON) {
				if (weStack_.size()>=1) { 
					dmrgWave_.we=weStack_.top();
					weStack_.pop();
					//std::cerr<<"CHANGED We taken from stack\n";
				} else {
					std::cerr<<"PUSHING STACK ERROR E\n";
					throw std::runtime_error("Environ Stack is empty\n");
				}
			}
			//std::cerr<<"PUSHING (POPPING) STACKSIZE="<<weStack_.size()<<" ";
			//std::cerr<<pSprime.block().size()<<"+"<<pEprime.block().size()<<"\n";
			if (counter_==0 && stage_==SHRINK_ENVIRON) {
// 				dmrgWave_.pEprime=pEprime;
// 				dmrgWave_.pSE=pSE;
// 				dmrgWave_.pSprime=pSprime;
// 				dmrgWave_.m=m;
// 				counter_++;
				//throw std::runtime_error("WFT::beforeWft(): Can't apply WFT\n");
				//return;
				if (weStack_.size()>=1) { 
					dmrgWave_.we=weStack_.top();
					//weStack_.pop();
					//std::cerr<<"CHANGED-COUNTER0 We taken from stack\n";
				} else {
					std::cerr<<"PUSHING-COUNTER0 STACK ERROR E\n";
					throw std::runtime_error("Environ Stack is empty\n");
				}
			}
			
			if (counter_==0 && stage_==SHRINK_SYSTEM) {
				//matrixIdentity(dmrgWave_.we,sizeOfOneSiteHilbertSpace_);
				//matrixIdentity(dmrgWave_.ws,dmrgWave_.ws.n_row());
// 				dmrgWave_.pEprime=pEprime;
// 				dmrgWave_.pSE=pSE;
// 				dmrgWave_.pSprime=pSprime;
// 				dmrgWave_.m=m;
// 				counter_++;
				//throw std::runtime_error("WFT::beforeWft(): Can't apply WFT\n");
				//return;
				
			}
		}
		
		template<typename SomeVectorType>
		void createVector(SomeVectorType& psiDest,const SomeVectorType& psiSrc,const BasisWithOperatorsType& pSprime,
				  const BasisWithOperatorsType& pEprime,const BasisType& pSE) const
		{
			if (useSu2Symmetry_)
				transformVectorSu2(psiDest,psiSrc,pSprime,pEprime,pSE);
			else
				transformVector(psiDest,psiSrc,pSprime,pEprime,pSE);
			
			RealType eps = 1e-4;
			if (std::norm(psiDest)<eps) {
				std::cerr<<"ATTENTION norm="<<std::norm(psiDest)<<" is too small\n";
				std::cerr<<"ATTENTION originalNorm="<<std::norm(psiSrc)<<"\n";
				//createRandomVector(psiDest);
				throw std::runtime_error(" createVector(): norm<1\n");
			}
			std::ostringstream msg;
			msg<<"I'm working hard!";
			progress_.printline(msg,std::cout);
		}
			
		void afterWft(const BasisWithOperatorsType& pSprime,
				  const BasisWithOperatorsType& pEprime,const BasisType& pSE) //,size_t m)
		{
			dmrgWave_.pEprime=pEprime;
			dmrgWave_.pSE=pSE;
			dmrgWave_.pSprime=pSprime;
			//dmrgWave_.m=m;
			firstCall_=false;
			counter_++;
		}
		
		template<typename SomeVectorType>
		void transformVector(SomeVectorType& psiDest,const SomeVectorType& psiSrc,const BasisWithOperatorsType& pSprime,
				      const BasisWithOperatorsType& pEprime,const BasisType& pSE) const
		{
			//std::cerr<<"counter="<<counter_<<"direction = "<<stage_<<"\n";
			if (stage_==SHRINK_SYSTEM) {
				if (firstCall_) throw std::runtime_error("WFT: This corner case is unimplmemented yet (sorry!)\n");
				else if (counter_==0) transformVector1bounce(psiDest,psiSrc,pSprime,pEprime,pSE);
				else transformVector1(psiDest,psiSrc,pSprime,pEprime,pSE);
			}
			if (stage_==SHRINK_ENVIRON) {
				if (firstCall_) transformVector2FromInfinite(psiDest,psiSrc,pSprime,pEprime,pSE);
 				else if (counter_==0) transformVector2bounce(psiDest,psiSrc,pSprime,pEprime,pSE);
				else transformVector2(psiDest,psiSrc,pSprime,pEprime,pSE);
			}
		}
		
		
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
			size_t nk = sizeOfOneSiteHilbertSpace_;
			size_t nip = pSE.permutationInverse().size()/pEprime.permutationInverse().size();
			size_t njp = pEprime.permutationInverse().size()/nk;
			//printDmrgWave();
			if (dmrgWave_.pSprime.permutationInverse().size()!=dmrgWave_.ws.n_row()) {
				printDmrgWave();
				throw std::runtime_error("transformVector1():"
						"SpermutationInverse.size()!=dmrgWave_.ws.n_row()\n");
			}
			if (njp!=dmrgWave_.we.n_col()) {
				printDmrgWave();
				std::cerr<<"nip="<<nip<<" njp="<<njp<<" nk="<<nk<<" dmrgWave_.we.n_col()="<<dmrgWave_.we.n_col()<<"\n";
				throw std::runtime_error("WaveFunctionTransformation::transformVector1():"
						"njp!=dmrgWave_.we.n_col()\n");
			}
			if (dmrgWave_.pSE.permutationInverse().size()!=psiSrc.size()) {
				printDmrgWave();
				std::cerr<<"SEpermutationInverse.size="<<dmrgWave_.pSE.permutationInverse().size();
				std::cerr<<" psiSrc.size="<<psiSrc.size()<<"\n";
				throw std::runtime_error("WaveFunctionTransformation::transformVector1():"
						" dmrgWave_.SEpermutationInverse.size()!=dmrgWave_.psi.size()\n");
			}
			
			size_t start = psiDest.offset(i0);
			size_t final = psiDest.effectiveSize(i0)+start;
			
			SparseMatrixType ws(dmrgWave_.ws);
			SparseMatrixType we(dmrgWave_.we);
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
			size_t nk = sizeOfOneSiteHilbertSpace_;
			size_t ni=dmrgWave_.ws.n_col();

			//int m = dmrgWave_.m;
			//size_t final = dmrgWave_.pSE.partition(m+1);
			//size_t start = dmrgWave_.pSE.partition(m);
			
			size_t nip = dmrgWave_.pSprime.permutationInverse().size()/nk;
			size_t alpha = dmrgWave_.pSprime.permutationInverse(ip+kp*nip);
			
			size_t offseta = 0; 
			size_t totala =dmrgWave_.ws.n_col();
			if (totala == dmrgWave_.pSprime.size()) {
				size_t eqn =  dmrgWave_.pSprime.qn(alpha,BasisType::BEFORE_TRANSFORM); 
				int ma = dmrgWave_.pSprime.partitionFromQn(eqn);
				if (ma>=0) {
					offseta = dmrgWave_.pSprime.partition(ma);
					totala = dmrgWave_.pSprime.partition(ma+1) - offseta;
				}
			}

			size_t offsetb = 0;
			size_t totalb = dmrgWave_.pEprime.permutationInverse().size();
			if  (dmrgWave_.pEprime.dmrgTransformed()) {
				size_t eqn = dmrgWave_.pEprime.qn(jp); 
				int mb =  dmrgWave_.pEprime.partitionFromQn(eqn,BasisType::BEFORE_TRANSFORM);
				if (mb>=0) {
					offsetb = dmrgWave_.pEprime.partition(mb,BasisType::BEFORE_TRANSFORM);
					totalb = dmrgWave_.pEprime.partition(mb+1,BasisType::BEFORE_TRANSFORM) - offsetb;
				}
			}
			
			SparseElementType sum=0;

			for (int k = ws.getRowPtr(alpha);k<ws.getRowPtr(alpha+1);k++) {
				size_t i = ws.getCol(k);
				for (int k2=weT.getRowPtr(jp);k2<weT.getRowPtr(jp+1);k2++) {
					size_t j = weT.getCol(k2);
					size_t x = dmrgWave_.pSE.permutationInverse(i+j*ni);
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
			size_t nk = sizeOfOneSiteHilbertSpace_;
			size_t nip = pSprime.permutationInverse().size()/nk;
			size_t nalpha = pSprime.permutationInverse().size();
			//printDmrgWave();
			if (dmrgWave_.pEprime.permutationInverse().size()!=dmrgWave_.we.n_row()) {
				printDmrgWave();
				throw std::runtime_error("transformVector2():"
						"PpermutationInverse.size()!=dmrgWave_.we.n_row()\n");
			}
			if (nip!=dmrgWave_.ws.n_col()) {
				printDmrgWave();
				throw std::runtime_error("WaveFunctionTransformation::transformVector2():"
						"nip!=dmrgWave_.ws.n_row()\n");
			}
			if (dmrgWave_.pSE.permutationInverse().size()!=psiSrc.size()) {
				printDmrgWave();
				std::cerr<<"SEpermutationInverse.size="<<dmrgWave_.pSE.permutationInverse().size();
				std::cerr<<" psiSrc.size="<<psiSrc.size()<<"\n";
				throw std::runtime_error("WaveFunctionTransformation::transformVector2():"
						" dmrgWave_.SEpermutationInverse.size()!=dmrgWave_.psi.size()\n");
			}

			size_t start = psiDest.offset(i0);
			size_t final = psiDest.effectiveSize(i0)+start;
			
			SparseMatrixType we(dmrgWave_.we);
			SparseMatrixType ws(dmrgWave_.ws);
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
			size_t nk = sizeOfOneSiteHilbertSpace_;
			size_t nalpha=dmrgWave_.pSprime.permutationInverse().size();
			
			size_t beta = dmrgWave_.pEprime.permutationInverse(kp+jp*nk);
			
			SparseElementType sum=0;
			
			for (int k=wsT.getRowPtr(ip);k<wsT.getRowPtr(ip+1);k++) {
				//SparseElementType sum2=0;
				size_t alpha = wsT.getCol(k);
				for (int k2=we.getRowPtr(beta);k2<we.getRowPtr(beta+1);k2++) {
					size_t j = we.getCol(k2);
					size_t x = dmrgWave_.pSE.permutationInverse(alpha+j*nalpha);
					sum += wsT.getValue(k)*we.getValue(k2)*psiSrc[x];
				}
				//sum += sum2;
			}
			return sum;
		}

		// SU(2) support below:
		template<typename SomeVectorType>
		void transformVectorSu2(SomeVectorType& psiDest,const SomeVectorType& psiSrc,const BasisWithOperatorsType& pSprime,
				      const BasisWithOperatorsType& pEprime,const BasisType& pSE) const
		{
			if (stage_==SHRINK_SYSTEM)
				transformVector1Su2(psiDest,psiSrc,pSprime,pEprime,pSE);
			if (stage_==SHRINK_ENVIRON)
				transformVector2Su2(psiDest,psiSrc,pSprime,pEprime,pSE);
		}
		
		template<typename SomeVectorType>
		void transformVector1Su2(SomeVectorType& psiDest,const SomeVectorType& psiSrc,const BasisWithOperatorsType& pSprime,
				      const BasisWithOperatorsType& pEprime,const BasisType& pSE) const
		{
			size_t nk = sizeOfOneSiteHilbertSpace_;
			size_t nip = pSE.getFactors().rank()/pEprime.getFactors().rank();
			size_t njp = pEprime.getFactors().rank()/nk;

			if ((size_t)dmrgWave_.pSprime.getFactors().rank()!=dmrgWave_.ws.n_row()) {
				printDmrgWave();
				throw std::runtime_error("transformVector1Su2(): getFactors.size()!=dmrgWave_.ws.n_row()\n");
			}
			if (njp!=dmrgWave_.we.n_col()) {
				printDmrgWave();
				std::cerr<<"nip="<<nip<<" njp="<<njp<<" nk="<<nk;
				std::cerr<<" dmrgWave_.we.n_col()="<<dmrgWave_.we.n_col()<<"\n";
				throw std::runtime_error("WaveFunctionTransformation::transformVector1Su2():"
						"njp!=dmrgWave_.we.n_col()\n");
			}
			if ((size_t)dmrgWave_.pSE.getFactors().rank()!=psiSrc.size()) {
				printDmrgWave();
				std::cerr<<"getFactors.size="<<dmrgWave_.pSE.permutationInverse().size();
				std::cerr<<" psiSrc.size="<<psiSrc.size()<<"\n";
				throw std::runtime_error("WaveFunctionTransformation::transformVector1Su2():"
						" dmrgWave_.getFactors.size()!=dmrgWave_.psi.size()\n");
			}

			for (size_t ii=0;ii<psiDest.sectors();ii++) {
				size_t i = psiDest.sector(ii);
				size_t start = psiDest.offset(i);
				size_t final = psiDest.effectiveSize(i)+start;
				transformVector1Su2(psiDest,psiSrc,pSprime,pEprime,pSE,start,final);
			}
		}
		
		template<typename SomeVectorType>
		void transformVector1Su2(SomeVectorType& psiDest,const SomeVectorType& psiSrc,
			const BasisWithOperatorsType& pSprime,
   			const BasisWithOperatorsType& pEprime,const BasisType& pSE,size_t start,size_t final) const
		{
			const FactorsType& factorsSE = pSE.getFactors();
			const FactorsType& factorsSEOld = dmrgWave_.pSE.getFactors();
			const FactorsType& factorsE = pEprime.getFactors();
			size_t nk = sizeOfOneSiteHilbertSpace_;
			size_t nip = pSE.getFactors().rank()/pEprime.getFactors().rank();
			
			FactorsType factorsInverseSE,
   				//factorsInverseSEOld,
   				factorsInverseE;
			transposeConjugate(factorsInverseSE,factorsSE);
			//transposeConjugate(factorsInverseSEOld,factorsSEOld);
			transposeConjugate(factorsInverseE,factorsE);
			
			SparseMatrixType ws(dmrgWave_.ws),we(dmrgWave_.we),weT;
			transposeConjugate(weT,we);
			
			for (size_t x=start;x<final;x++) {
				psiDest[x] = 0;
				for (int kI = factorsInverseSE.getRowPtr(x);kI < factorsInverseSE.getRowPtr(x+1);kI++) {
					size_t ip,beta;
					utils::getCoordinates(ip,beta,(size_t)factorsInverseSE.getCol(kI),nip);
					for (int k2I = factorsInverseE.getRowPtr(beta);k2I < factorsInverseE.getRowPtr(beta+1);k2I++) {
						size_t kp,jp;
						utils::getCoordinates(kp,jp,(size_t)factorsInverseE.getCol(k2I),nk);
						psiDest[x] += createVectorAux1bSu2(psiSrc,ip,kp,jp,factorsSEOld,ws,weT)*
								factorsInverseSE.getValue(kI)*factorsInverseE.getValue(k2I);
					}
				}
			}
		}
		
		template<typename SomeVectorType>
		SparseElementType createVectorAux1bSu2(const SomeVectorType& psiSrc,
					size_t ip,size_t kp,size_t jp, const FactorsType& factorsSE,
					      const SparseMatrixType& ws,const SparseMatrixType& weT) const
		{
			size_t nk = sizeOfOneSiteHilbertSpace_;
			size_t ni=dmrgWave_.ws.n_col();
			const FactorsType& factorsS = dmrgWave_.pSprime.getFactors();
			SparseElementType sum=0;
			size_t nip = dmrgWave_.pSprime.permutationInverse().size()/nk;
			
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
		void transformVector2Su2(SomeVectorType& psiDest,const SomeVectorType& psiSrc,const BasisWithOperatorsType& pSprime,
				      const BasisWithOperatorsType& pEprime,const BasisType& pSE) const
		{
			size_t nk = sizeOfOneSiteHilbertSpace_;
			size_t nip = pSprime.getFactors().rank()/nk;
			
			if (dmrgWave_.ws.n_row()!=dmrgWave_.pSprime.permutationInverse().size()) throw std::runtime_error("Error!!");
			if (dmrgWave_.we.n_col()!=dmrgWave_.pEprime.size()) throw std::runtime_error("Error\n");

			if ((size_t)dmrgWave_.pEprime.getFactors().rank()!=dmrgWave_.we.n_row()) {
				printDmrgWave();
				throw std::runtime_error("transformVector2Su2():"
						"PpermutationInverse.size()!=dmrgWave_.we.n_row()\n");
			}
			if (nip!=dmrgWave_.ws.n_col()) {
				printDmrgWave();
				throw std::runtime_error("WaveFunctionTransformation::transformVector2Su2():"
						"nip!=dmrgWave_.ws.n_row()\n");
			}
			if (dmrgWave_.pSE.permutationInverse().size()!=psiSrc.size()) {
				printDmrgWave();
				std::cerr<<"SEpermutationInverse.size="<<dmrgWave_.pSE.permutationInverse().size();
				std::cerr<<" psiSrc.size="<<psiSrc.size()<<"\n";
				throw std::runtime_error("WaveFunctionTransformation::transformVector2Su2():"
						" dmrgWave_.SEpermutationInverse.size()!=dmrgWave_.psi.size()\n");
			}

			for (size_t ii=0;ii<psiDest.sectors();ii++) {
				size_t i = psiDest.sector(ii);
				size_t start = psiDest.offset(i);
				size_t final = psiDest.effectiveSize(i)+start;
				transformVector2Su2(psiDest,psiSrc,pSprime,
				      pEprime,pSE,start,final);
			}
		}
		
		template<typename SomeVectorType>
		void transformVector2Su2(SomeVectorType& psiDest,const SomeVectorType& psiSrc,const BasisWithOperatorsType& pSprime,
				      const BasisWithOperatorsType& pEprime,const BasisType& pSE,size_t start,size_t final) const
		{
			size_t nk = sizeOfOneSiteHilbertSpace_;
			size_t nip = pSprime.getFactors().rank()/nk;
			size_t nalpha = pSprime.getFactors().rank();
			
			const FactorsType& factorsSE = pSE.getFactors();
			//const SparseMatrixType& factorsSEOld = dmrgWave_.pSE.getFactors();
			const FactorsType& factorsS = pSprime.getFactors();
			FactorsType factorsInverseSE,
   				//factorsInverseSEOld,
   					factorsInverseS;
			transposeConjugate(factorsInverseSE,factorsSE);
			//transposeConjugate(factorsInverseSEOld,factorsSEOld);
			transposeConjugate(factorsInverseS,factorsS);
			SparseMatrixType ws(dmrgWave_.ws),we(dmrgWave_.we),wsT;
			transposeConjugate(wsT,ws);
			
			for (size_t x=start;x<final;x++) {
				psiDest[x] = 0;
				size_t xx = x; // pSE.permutationInverse(x);
				for (int kI=factorsInverseSE.getRowPtr(xx);kI<factorsInverseSE.getRowPtr(xx+1);kI++) {
					size_t alpha,jp;
					utils::getCoordinates(alpha,jp,(size_t)factorsInverseSE.getCol(kI),nalpha);
					size_t alphax =  alpha; //pSprime.permutationInverse(alpha);
					for (int k2I=factorsInverseS.getRowPtr(alphax);k2I<factorsInverseS.getRowPtr(alphax+1);k2I++) {
						size_t ip,kp;
						utils::getCoordinates(ip,kp,(size_t)factorsInverseS.getCol(k2I),nip);
						psiDest[x] += fastAux2bSu2(psiSrc,ip,kp,jp,wsT,we)* //factorsInverseSEOld)*
								factorsInverseSE.getValue(kI)*factorsInverseS.getValue(k2I);
					}
				}
			}
		}
		
		template<typename SomeVectorType>
		SparseElementType fastAux2bSu2(const SomeVectorType& psiSrc,size_t ip,size_t kp,size_t jp,
					const SparseMatrixType& wsT,const SparseMatrixType& we) const
		{
			size_t nk = sizeOfOneSiteHilbertSpace_;
			size_t nalpha=dmrgWave_.pSprime.getFactors().rank();
			SparseElementType sum=0;
			const FactorsType& factorsE = dmrgWave_.pEprime.getFactors();
			const FactorsType& factorsSE = dmrgWave_.pSE.getFactors();
			//int m = dmrgWave_.m;
			//size_t final = dmrgWave_.pSE.partition(m+1);
			//size_t start = dmrgWave_.pSE.partition(m);
			size_t eqn =  dmrgWave_.pSprime.qn(ip);
			
			int ma =  dmrgWave_.pSprime.partitionFromQn(eqn,BasisType::BEFORE_TRANSFORM);
			if (ma<0) throw std::runtime_error("ma<0\n");
			
			size_t totala =  dmrgWave_.ws.n_row();
			totala = dmrgWave_.pSprime.partition(ma+1,BasisType::BEFORE_TRANSFORM)-
					dmrgWave_.pSprime.partition(ma,BasisType::BEFORE_TRANSFORM);
			//size_t offseta = dmrgWave_.pSprime.partition(ma,DmrgBasisType::BEFORE_TRANSFORM);
			
			size_t kpjp = kp+jp*nk;
			size_t kpjpx = dmrgWave_.pEprime.permutationInverse(kpjp);
			
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

		void printDmrgWave() const
		{
			std::cerr<<dmrgWave_;
			std::cerr<<"wsStack="<<wsStack_.size()<<"\n";
			std::cerr<<"weStack="<<weStack_.size()<<"\n";
			std::cerr<<"counter="<<counter_<<"\n"; 
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
			size_t nk = sizeOfOneSiteHilbertSpace_;
			size_t nip = pSprime.permutationInverse().size()/nk;
			size_t nalpha = pSprime.permutationInverse().size();
			
			std::ostringstream msg;
			msg<<" We're moving to the finite loop, bumpy ride ahead!";
			progress_.printline(msg,std::cout);
			
			if (dmrgWave_.pEprime.permutationInverse().size()!=dmrgWave_.we.n_row()) {
				printDmrgWave();
				throw std::runtime_error("transformVector2():"
						"PpermutationInverse.size()!=dmrgWave_.we.n_row()\n");
			}
			if (nip!=dmrgWave_.ws.n_col()) {
				printDmrgWave();
				throw std::runtime_error("WaveFunctionTransformation::transformVector2():"
						"nip!=dmrgWave_.ws.n_row()\n");
			}
			if (dmrgWave_.pSE.permutationInverse().size()!=psiSrc.size()) {
				printDmrgWave();
				std::cerr<<"SEpermutationInverse.size="<<dmrgWave_.pSE.permutationInverse().size();
				std::cerr<<" psiSrc.size="<<psiSrc.size()<<"\n";
				throw std::runtime_error("WaveFunctionTransformation::transformVector2():"
						" dmrgWave_.SEpermutationInverse.size()!=dmrgWave_.psi.size()\n");
			}

			size_t start = psiDest.offset(i0);
			size_t final = psiDest.effectiveSize(i0)+start;
			
			SparseMatrixType we(dmrgWave_.we);
			SparseMatrixType ws(dmrgWave_.ws);
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
			size_t nk = sizeOfOneSiteHilbertSpace_;
			size_t nalpha=dmrgWave_.pSprime.permutationInverse().size();
			
			
			
			SparseElementType sum=0;
			
			for (int k=wsT.getRowPtr(is);k<wsT.getRowPtr(is+1);k++) {
				//SparseElementType sum2=0;
				size_t ip = wsT.getCol(k);
				SparseElementType sum2 = 0;
				//for (int k2=we.getRowPtr(jen);k2<we.getRowPtr(jen+1);k2++) {
					size_t jpr = jen; //we.getCol(k2);
					size_t jp = dmrgWave_.pEprime.permutationInverse(jpl + jpr*nk);
					size_t y = dmrgWave_.pSE.permutationInverse(ip + jp*nalpha);
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
			size_t nk = sizeOfOneSiteHilbertSpace_;
			size_t nip = pSE.permutationInverse().size()/pEprime.permutationInverse().size();
			//size_t njp = pEprime.permutationInverse().size()/nk;
			//printDmrgWave();
			std::ostringstream msg;
			msg<<" We're bouncing on the right, so buckle up!";
			progress_.printline(msg,std::cout);
			
			if (dmrgWave_.pSE.permutationInverse().size()!=psiSrc.size()) {
				printDmrgWave();
				std::cerr<<"SEpermutationInverse.size="<<dmrgWave_.pSE.permutationInverse().size();
				std::cerr<<" psiSrc.size="<<psiSrc.size()<<"\n";
				throw std::runtime_error("WaveFunctionTransformation::transformVector1():"
						" dmrgWave_.SEpermutationInverse.size()!=dmrgWave_.psi.size()\n");
			}
			
			size_t start = psiDest.offset(i0);
			size_t final = psiDest.effectiveSize(i0)+start;
			
			/* SparseMatrixType ws(dmrgWave_.ws);
			SparseMatrixType we(dmrgWave_.we);
			SparseMatrixType weT;
			transposeConjugate(weT,we);
			*/
			//std::cerr<<"ALTERNATIVE\n";
			size_t nalpha=dmrgWave_.pSprime.permutationInverse().size();
			for (size_t x=start;x<final;x++) {
				size_t ip,beta,kp,jp;
				utils::getCoordinates(ip,beta,(size_t)pSE.permutation(x),nip);
				utils::getCoordinates(kp,jp,(size_t)pEprime.permutation(beta),nk);
				size_t ipkp = dmrgWave_.pSprime.permutationInverse(ip + kp*nip);
				size_t y = dmrgWave_.pSE.permutationInverse(ipkp + jp*nalpha);
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
			size_t nk = sizeOfOneSiteHilbertSpace_;
			size_t nip = pSprime.permutationInverse().size()/nk;
			size_t nalpha = pSprime.permutationInverse().size();
			//printDmrgWave();
			
			std::ostringstream msg;
			msg<<" We're bouncing on the left, so buckle up!";
			progress_.printline(msg,std::cout);
			
			if (dmrgWave_.pSE.permutationInverse().size()!=psiSrc.size()) {
				printDmrgWave();
				std::cerr<<"SEpermutationInverse.size="<<dmrgWave_.pSE.permutationInverse().size();
				std::cerr<<" psiSrc.size="<<psiSrc.size()<<"\n";
				throw std::runtime_error("WaveFunctionTransformation::transformVector2():"
						" dmrgWave_.SEpermutationInverse.size()!=dmrgWave_.psi.size()\n");
			}

			size_t start = psiDest.offset(i0);
			size_t final = psiDest.effectiveSize(i0)+start;
			
			
			for (size_t x=start;x<final;x++) {
				size_t ip,alpha,kp,jp;
				utils::getCoordinates(alpha,jp,(size_t)pSE.permutation(x),nalpha);
				utils::getCoordinates(ip,kp,(size_t)pSprime.permutation(alpha),nip);
				size_t kpjp = dmrgWave_.pEprime.permutationInverse(kp + jp*nk);
				
				size_t y = dmrgWave_.pSE.permutationInverse(ip + kpjp*nip);
				psiDest[x]=psiSrc[y];
			}
			
		}

	}; // class WaveFunctionTransformation
} // namespace Dmrg

/*@}*/
#endif
