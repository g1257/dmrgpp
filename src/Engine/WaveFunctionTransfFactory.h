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
 *  This class implements the wave function transformation factory,
 *  see PRL 77, 3633 (1996)
 *
 */

#ifndef WFT_FACTORY_H
#define WFT_FACTORY_H
 
#include "ProgressIndicator.h"
#include "WaveFunctionTransfLocal.h"
#include "WaveFunctionTransfSu2.h"
#include "DmrgWaveStruct.h"
#include "IoSimple.h"

namespace Dmrg {
	


	template<typename LeftRightSuperType,typename VectorWithOffsetType>
	class WaveFunctionTransfFactory {
		typedef PsimagLite::IoSimple IoType;
		public:
		enum {DO_NOT_RESET_COUNTER,RESET_COUNTER};

		typedef typename LeftRightSuperType::BasisWithOperatorsType
			BasisWithOperatorsType;
		typedef typename BasisWithOperatorsType::SparseMatrixType SparseMatrixType;
		typedef typename BasisWithOperatorsType::BasisType BasisType;
		typedef typename SparseMatrixType::value_type SparseElementType;
		typedef std::vector<SparseElementType> VectorType;
		typedef typename BasisWithOperatorsType::RealType RealType;
		typedef typename BasisType::FactorsType FactorsType;
		typedef DmrgWaveStruct<LeftRightSuperType> DmrgWaveStructType;

		typedef WaveFunctionTransfBase<DmrgWaveStructType,VectorWithOffsetType>
					WaveFunctionTransfBaseType;
		typedef WaveFunctionTransfLocal<DmrgWaveStructType,VectorWithOffsetType>
			WaveFunctionTransfLocalType;
		typedef WaveFunctionTransfSu2<DmrgWaveStructType,VectorWithOffsetType>
					WaveFunctionTransfSu2Type;

		static const size_t INFINITE = ProgramGlobals::INFINITE;
		static const size_t EXPAND_SYSTEM = ProgramGlobals::EXPAND_SYSTEM;
		static const size_t EXPAND_ENVIRON = ProgramGlobals::EXPAND_ENVIRON;
		
		template<typename SomeParametersType>
		WaveFunctionTransfFactory(SomeParametersType& parameters,size_t nk)
		: hilbertSpaceOneSite_(nk),
		  isEnabled_(!(parameters.options.find("nowft")!=std::string::npos)),
		  stage_(INFINITE),
		  counter_(0),
		  firstCall_(true),
		  progress_("WaveFunctionTransf",0),
		  filenameIn_(parameters.checkpoint.filename),
		  filenameOut_(parameters.filename),
		  WFT_STRING("Wft"),
		  wftImpl_(0)
		{
			if (!isEnabled_) return;
			if (parameters.options.find("checkpoint")!=std::string::npos)
				load();
			if (BasisType::useSu2Symmetry()) {
				wftImpl_=new WaveFunctionTransfSu2Type(hilbertSpaceOneSite_,
						stage_,firstCall_,counter_,dmrgWaveStruct_);
			} else {
				wftImpl_=new WaveFunctionTransfLocalType(hilbertSpaceOneSite_,
						stage_,firstCall_,counter_,dmrgWaveStruct_);
			}
		}

		~WaveFunctionTransfFactory()
		{
			if (!isEnabled_) return;
			save();
			delete wftImpl_;
		}

		void setStage(size_t stage)
		{
			if (stage == stage_) return;
			stage_=stage;
			counter_=0;
		}
		
		void triggerOn(const LeftRightSuperType& lrs)
		{
			bool allow=false;
			switch (stage_) {
				case INFINITE:
					allow=false;
					break;
				case EXPAND_SYSTEM:
					allow=true;
					
				case EXPAND_ENVIRON:
					allow=true;
			}
			// FIXME: Must check the below change when using SU(2)!!
			//if (m<0) allow = false; // isEnabled_=false;
			
			if (!isEnabled_ || !allow) return;
			//try {
				beforeWft(lrs);
			//} catch (std::exception& e) {
			//	doNextOne_=false;
			//}
			std::ostringstream msg;
			msg<<"Window open, ready to transform vectors";
			progress_.printline(msg,std::cout);
		}
		
		// FIXME: change name to transformVector
		template<typename SomeVectorType,typename SomeVectorType2>
		void setInitialVector(	
					SomeVectorType& dest,
					const SomeVectorType2& src,
					const LeftRightSuperType& lrs) const
		{
			bool allow=false;
			switch (stage_) {
				case INFINITE:
					allow=false;
					break;
				case EXPAND_SYSTEM:
					allow=true;
					
				case EXPAND_ENVIRON:
					allow=true;
			}
			// FIXME: Must check the below change when using SU(2)!!
			//if (m<0) allow = false; // isEnabled_=false;
			
			if (isEnabled_ && allow) {
				RealType eps = 1e-6;
				if (std::norm(src)<eps)
					throw std::runtime_error("src's norm is zero\n");
				createVector(dest,src,lrs);
			} else {
				createRandomVector(dest);
			}
		}
		
		void triggerOff(const LeftRightSuperType& lrs)
		{
			bool allow=false;
			switch (stage_) {
				case INFINITE:
					allow=false;
					break;
				case EXPAND_SYSTEM:
					allow=true;
					
				case EXPAND_ENVIRON:
					allow=true;
			}
			// FIXME: Must check the below change when using SU(2)!!
			//if (m<0) allow = false; // isEnabled_=false;
			
			if (!isEnabled_ || !allow) return;
			afterWft(lrs);
			//doNextOne_=true;
			std::ostringstream msg;
			msg<<"Window closed, no more transformations, please";
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
				myRandomT(tmp);
				y[i]=tmp;
				atmp += myProductT(y[i],y[i]);
			}
			atmp = 1.0 / sqrt (atmp);
			for (size_t i=offset;i<final;i++) y[i] *= atmp;

		}

		template<typename SomeMatrixType>
		void push(
			const SomeMatrixType& transform,
			size_t direction,
			LeftRightSuperType& lrs)
		{
			if (!isEnabled_) return;
			
			switch (stage_) {
				case INFINITE:
					if (direction==EXPAND_SYSTEM) {
						wsStack_.push(transform);
						dmrgWaveStruct_.ws=transform;
					} else {
						weStack_.push(transform);
						dmrgWaveStruct_.we=transform;
						//std::cerr<<"CHANGED dmrgWaveStruct_.we to transform\n";
						//std::cerr<<"PUSHING "<<transform.n_row()<<"x"<<transform.n_col()<<"\n";
					}
					break;
				case EXPAND_ENVIRON:
					if (direction!=EXPAND_ENVIRON) throw std::logic_error("EXPAND_ENVIRON but option==0\n");
					dmrgWaveStruct_.we=transform;
					dmrgWaveStruct_.ws=transform;
					//vectorConvert(dmrgWaveStruct_.psi,psi);
					weStack_.push(transform);
					//std::cerr<<"PUSHING (POPPING) We "<<weStack_.size()<<"\n";
					break;
				case EXPAND_SYSTEM:
					if (direction!=EXPAND_SYSTEM) throw std::logic_error("EXPAND_SYSTEM but option==1\n");
					dmrgWaveStruct_.ws=transform;
					dmrgWaveStruct_.we=transform;
					//vectorConvert(dmrgWaveStruct_.psi,psi);
					wsStack_.push(transform);
					break;
			}

			dmrgWaveStruct_.lrs=lrs;
//			if (direction==EXPAND_SYSTEM) { // transforming the system
//				dmrgWaveStruct_.pEprime=pBasisSummed;
//				dmrgWaveStruct_.pSprime=pBasis;
//			} else {
//				dmrgWaveStruct_.pSprime=pBasisSummed;
//				dmrgWaveStruct_.pEprime=pBasis;
//			}
			std::ostringstream msg;
			msg<<"OK, pushing option="<<direction<<" and stage="<<stage_;
			progress_.printline(msg,std::cout);
		}

//		void disable() { isEnabled_=false; }
//
		bool isEnabled() const { return isEnabled_; }
		
	private:
		
		void beforeWft(const LeftRightSuperType& lrs)
		{
			if (stage_==EXPAND_ENVIRON) {
				if (wsStack_.size()>=1) {
					dmrgWaveStruct_.ws=wsStack_.top();
					wsStack_.pop();
				} else {
					//std::cerr<<"PUSHING STACK ERROR S\n";
					throw std::runtime_error("System Stack is empty\n");
				}
			}
			
			if (stage_==EXPAND_SYSTEM) {
				if (weStack_.size()>=1) { 
					dmrgWaveStruct_.we=weStack_.top();
					weStack_.pop();
					//std::cerr<<"CHANGED We taken from stack\n";
				} else {
					//std::cerr<<"PUSHING STACK ERROR E\n";
					throw std::runtime_error("Environ Stack is empty\n");
				}
			}
			//std::cerr<<"PUSHING (POPPING) STACKSIZE="<<weStack_.size()<<" ";
			//std::cerr<<pSprime.block().size()<<"+"<<pEprime.block().size()<<"\n";
			if (counter_==0 && stage_==EXPAND_SYSTEM) {
// 				dmrgWaveStruct_.pEprime=pEprime;
// 				dmrgWaveStruct_.pSE=pSE;
// 				dmrgWaveStruct_.pSprime=pSprime;
// 				dmrgWaveStruct_.m=m;
// 				counter_++;
				//throw std::runtime_error("WFT::beforeWft(): Can't apply WFT\n");
				//return;
				if (weStack_.size()>=1) { 
					dmrgWaveStruct_.we=weStack_.top();
					//weStack_.pop();
					//std::cerr<<"CHANGED-COUNTER0 We taken from stack\n";
				} else {
					std::cerr<<"PUSHING-COUNTER0 STACK ERROR E\n";
					throw std::runtime_error("Environ Stack is empty\n");
				}
			}
			
			if (counter_==0 && stage_==EXPAND_ENVIRON) {
				//matrixIdentity(dmrgWaveStruct_.we,hilbertSpaceOneSite_);
				//matrixIdentity(dmrgWaveStruct_.ws,dmrgWaveStruct_.ws.n_row());
// 				dmrgWaveStruct_.pEprime=pEprime;
// 				dmrgWaveStruct_.pSE=pSE;
// 				dmrgWaveStruct_.pSprime=pSprime;
// 				dmrgWaveStruct_.m=m;
// 				counter_++;
				//throw std::runtime_error("WFT::beforeWft(): Can't apply WFT\n");
				//return;
				if (wsStack_.size()>=1) {
					dmrgWaveStruct_.ws=wsStack_.top();
					//weStack_.pop();
					//std::cerr<<"CHANGED-COUNTER0 We taken from stack\n";
				} else {
					std::cerr<<"PUSHING-COUNTER0 STACK ERROR E\n";
					throw std::runtime_error("System Stack is empty\n");
				}
			}
		}
		
		void createVector(
				VectorWithOffsetType& psiDest,
				const VectorWithOffsetType& psiSrc,
				const LeftRightSuperType& lrs) const
		{
			//try {
				wftImpl_->transformVector(psiDest,psiSrc,lrs);
			//} catch (std::exception& e) {
			//	printDmrgWave();
			//	throw e;
			//}
			std::ostringstream msg;
			msg<<"Transformation completed";
			progress_.printline(msg,std::cout);
		}
			
		void afterWft(const LeftRightSuperType& lrs)
		{
			dmrgWaveStruct_.lrs=lrs;
			firstCall_=false;
			counter_++;
		}

		void printDmrgWave() const
		{
			PsimagLite::IoSimple::Out io(std::cerr);
			dmrgWaveStruct_.save(io);
			std::cerr<<"wsStack="<<wsStack_.size()<<"\n";
			std::cerr<<"weStack="<<weStack_.size()<<"\n";
			std::cerr<<"counter="<<counter_<<"\n";
		}

		void save() const
		{
			if (!isEnabled_) throw std::runtime_error(
					"WFT::save(...) called but wft is disabled\n");

			typename IoType::Out io(WFT_STRING + filenameOut_,0);
			std::string s="isEnabled="+ttos(isEnabled_);
			io.printline(s);
			s="stage="+ttos(stage_);
			io.printline(s);
			s="counter="+ttos(counter_);
			io.printline(s);
			io.printline("dmrgWaveStruct");
			dmrgWaveStruct_.save(io);
			io.printMatrix(wsStack_,"wsStack");
			io.printMatrix(weStack_,"weStack");
		}

		void load()
		{
			if (!isEnabled_) throw std::runtime_error(
					"WFT::load(...) called but wft is disabled\n");

			typename IoType::In io(WFT_STRING + filenameIn_);
			io.readline(isEnabled_,"isEnabled=");
			io.readline(stage_,"stage=");
			io.readline(counter_,"counter=");
			firstCall_=false;
			io.advance("dmrgWaveStruct");
			dmrgWaveStruct_.load(io);
			io.readMatrix(wsStack_,"wsStack");
			io.readMatrix(weStack_,"weStack");
		}

		template<typename Field>
		void myRandomT(std::complex<Field> &value)
		{
			value = std::complex<Field>(drand48 () - 0.5, drand48 () - 0.5);
		}

		template<typename Field>
		void myRandomT(Field &value)
		{
			value = drand48 () - 0.5;
		}

		template<typename Field>
		Field myProductT(Field const &value1,Field const &value2)
		{
			return std::real(value1*std::conj(value2));
		}

		template<typename Field>
		Field myProductT(std::complex<Field> const &value1,
				std::complex<Field> const &value2)
		{
			return real(value1*conj(value2));
		}

		size_t hilbertSpaceOneSite_;
		bool isEnabled_;
		size_t stage_;
		size_t counter_;
		bool firstCall_;
		PsimagLite::ProgressIndicator progress_;
		std::string filenameIn_,filenameOut_;
		const std::string WFT_STRING;
		DmrgWaveStructType dmrgWaveStruct_;
		std::stack<PsimagLite::Matrix<SparseElementType> > wsStack_,weStack_;
		WaveFunctionTransfBaseType* wftImpl_;
	}; // class WaveFunctionTransformation
} // namespace Dmrg

/*@}*/
#endif
