/*
Copyright (c) 2009-2012, UT-Battelle, LLC
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

/*! \file WaveFunctionTransformation.h
 *
 *  This class implements the wave function transformation factory,
 *  see PRL 77, 3633 (1996)
 *
 */

#ifndef WFT_FACTORY_H
#define WFT_FACTORY_H
#include "Utils.h"
#include "ProgressIndicator.h"
#include "WaveFunctionTransfLocal.h"
#include "WaveFunctionTransfSu2.h"
#include "DmrgWaveStruct.h"
#include "IoSimple.h"
#include "Random48.h"

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
		WaveFunctionTransfFactory(SomeParametersType& params)
		: isEnabled_(!(params.options.find("nowft")!=std::string::npos)),
		  stage_(INFINITE),
		  counter_(0),
		  firstCall_(true),
		  progress_("WaveFunctionTransf",0),
		  filenameIn_(params.checkpoint.filename),
		  filenameOut_(params.filename),
		  WFT_STRING("Wft"),
		  wftImpl_(0),
		  rng_(3433117)
		{
			if (!isEnabled_) return;
			if (params.options.find("checkpoint")!=std::string::npos || params.options.find("restart")!=std::string::npos)
				load();
			if (BasisType::useSu2Symmetry()) {
				wftImpl_=new WaveFunctionTransfSu2Type(stage_,firstCall_,counter_,dmrgWaveStruct_);
			} else {
				wftImpl_=new WaveFunctionTransfLocalType(stage_,firstCall_,counter_,dmrgWaveStruct_);
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
			beforeWft(lrs);
			std::ostringstream msg;
			msg<<"Window open, ready to transform vectors";
			progress_.printline(msg,std::cout);
		}

		// FIXME: change name to transformVector
		template<typename SomeVectorType,typename SomeVectorType2>
		void setInitialVector(	
					SomeVectorType& dest,
					const SomeVectorType2& src,
					const LeftRightSuperType& lrs,
					size_t nk) const
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
#ifndef NDEBUG
				RealType eps = 1e-12;
				RealType x =std::norm(src);
				bool b = (x<eps);
				if (b) std::cerr<<"norm="<<x<<"\n";
				assert(!b);
#endif
				assert(nk>0);
				createVector(dest,src,lrs,nk);
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
				atmp += std::real(y[i]*std::conj(y[i]));
			}
			atmp = 1.0 / sqrt (atmp);
			for (size_t i=offset;i<final;i++) y[i] *= atmp;

		}

		void push(const SparseMatrixType& transform,
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
					}
					break;
				case EXPAND_ENVIRON:
					if (direction!=EXPAND_ENVIRON) throw std::logic_error("EXPAND_ENVIRON but option==0\n");
					dmrgWaveStruct_.we=transform;
					dmrgWaveStruct_.ws=transform;
					weStack_.push(transform);
					break;
				case EXPAND_SYSTEM:
					if (direction!=EXPAND_SYSTEM) throw std::logic_error("EXPAND_SYSTEM but option==1\n");
					dmrgWaveStruct_.ws=transform;
					dmrgWaveStruct_.we=transform;
					wsStack_.push(transform);
					break;
			}

			dmrgWaveStruct_.lrs=lrs;
			std::ostringstream msg;
			msg<<"OK, pushing option="<<direction<<" and stage="<<stage_;
			progress_.printline(msg,std::cout);
		}

		const SparseMatrixType& transform(size_t what) const
		{
			return (what==ProgramGlobals::SYSTEM) ? dmrgWaveStruct_.ws : dmrgWaveStruct_.we;
		}

		bool isEnabled() const { return isEnabled_; }
		
	private:
		
		void myRandomT(std::complex<RealType> &value) const
		{
			value = std::complex<RealType>(rng_() - 0.5, rng_() - 0.5);
		}

		void myRandomT(RealType &value) const
		{
			value = rng_() - 0.5;
		}

		void beforeWft(const LeftRightSuperType& lrs)
		{
			if (stage_==EXPAND_ENVIRON) {
				if (wsStack_.size()>=1) {
					dmrgWaveStruct_.ws=wsStack_.top();
					wsStack_.pop();
				} else {
					throw std::runtime_error("System Stack is empty\n");
				}
			}
			
			if (stage_==EXPAND_SYSTEM) {
				if (weStack_.size()>=1) { 
					dmrgWaveStruct_.we=weStack_.top();
					weStack_.pop();
				} else {
					throw std::runtime_error("Environ Stack is empty\n");
				}
			}
			if (counter_==0 && stage_==EXPAND_SYSTEM) {
				if (weStack_.size()>=1) { 
					dmrgWaveStruct_.we=weStack_.top();
				} 
			}
			
			if (counter_==0 && stage_==EXPAND_ENVIRON) {
				if (wsStack_.size()>=1) {
					dmrgWaveStruct_.ws=wsStack_.top();
				} 
			}
		}
		
		void createVector(
				VectorWithOffsetType& psiDest,
				const VectorWithOffsetType& psiSrc,
				const LeftRightSuperType& lrs,
				size_t nk) const
		{
			wftImpl_->transformVector(psiDest,psiSrc,lrs,nk);

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
			io.print("isEnabled=",isEnabled_);

			io.print("stage=",stage_);

			io.print("counter=",counter_);

			io.print("dmrgWaveStruct");
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

		bool isEnabled_;
		size_t stage_;
		size_t counter_;
		bool firstCall_;
		PsimagLite::ProgressIndicator progress_;
		std::string filenameIn_,filenameOut_;
		const std::string WFT_STRING;
		DmrgWaveStructType dmrgWaveStruct_;
		std::stack<SparseMatrixType> wsStack_,weStack_;
		WaveFunctionTransfBaseType* wftImpl_;
		PsimagLite::Random48<RealType> rng_;
	}; // class WaveFunctionTransformation
} // namespace Dmrg

/*@}*/
#endif
