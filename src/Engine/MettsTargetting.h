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

#ifndef METTS_TARGETTING_H
#define METTS_TARGETTING_H
#include <iostream>
#include "ProgressIndicator.h"
#include "BLAS.h"
#include "ApplyOperatorLocal.h"
#include "MettsSerializer.h"
#include "MettsParams.h"
#include "MettsStochastics.h"
#include <cassert>
#include "MettsCollapse.h"
#include "VectorWithOffset.h"
#include "ParametersForSolver.h"
#include "RandomForTests.h"
#include "TimeSerializer.h"
#include "TimeVectorsKrylov.h"
#include "TimeVectorsRungeKutta.h"
#include "TimeVectorsSuzukiTrotter.h"

namespace Dmrg {
	template<
			template<typename,typename,typename> class LanczosSolverTemplate,
			template<typename,typename> class InternalProductTemplate,
	 		template<typename,typename> class WaveFunctionTransfTemplate,
    			typename ModelType_,
    			typename IoType_,
       			template<typename> class VectorWithOffsetTemplate>
	class MettsTargetting  {

			struct MettsPrev {
				MettsPrev() : fixed(0),permutationInverse(0) { }
				SizeType fixed;
				typename PsimagLite::Vector<SizeType>::Type permutationInverse;
			};

		public:
			typedef ModelType_ ModelType;
			typedef IoType_ IoType;
			typedef typename ModelType::RealType RealType;
			typedef InternalProductTemplate<RealType,ModelType>
			   InternalProductType;
			typedef typename ModelType::OperatorsType OperatorsType;
			typedef typename ModelType::ModelHelperType ModelHelperType;
			typedef typename ModelHelperType::LeftRightSuperType
			    LeftRightSuperType;
			typedef typename LeftRightSuperType::BasisWithOperatorsType
			    BasisWithOperatorsType;
			typedef typename BasisWithOperatorsType::BasisType BasisType;
			typedef typename PsimagLite::Vector<RealType>::Type VectorType;
			typedef PsimagLite::ParametersForSolver<RealType> ParametersForSolverType;
			typedef LanczosSolverTemplate<ParametersForSolverType,InternalProductType,VectorType>
			    LanczosSolverType;
			typedef typename LanczosSolverType::TridiagonalMatrixType TridiagonalMatrixType;
			typedef typename BasisWithOperatorsType::OperatorType OperatorType;
			typedef MettsParams<ModelType> TargettingParamsType;
			typedef typename BasisType::BlockType BlockType;
			typedef VectorWithOffsetTemplate<RealType> VectorWithOffsetType;
			typedef WaveFunctionTransfTemplate<LeftRightSuperType,VectorWithOffsetType> WaveFunctionTransfType;
			typedef VectorType TargetVectorType;
			typedef PsimagLite::Matrix<RealType> MatrixType;
			typedef BlockMatrix<RealType,MatrixType> BlockMatrixType;
			typedef ApplyOperatorLocal<LeftRightSuperType,VectorWithOffsetType,TargetVectorType> ApplyOperatorType;
			typedef MettsSerializer<RealType,VectorWithOffsetType> MettsSerializerType;
			typedef typename PsimagLite::RandomForTests<RealType> RngType;
			typedef MettsStochastics<ModelType,RngType> MettsStochasticsType;
			typedef typename MettsStochasticsType::PairType PairType;
			typedef MettsCollapse<VectorWithOffsetType,MettsStochasticsType,TargettingParamsType> MettsCollapseType;
			typedef typename MettsCollapseType::PackIndicesType PackIndicesType;
			typedef TimeSerializer<RealType,VectorWithOffsetType> TimeSerializerType;
			typedef TimeVectorsBase<TargettingParamsType,ModelType,WaveFunctionTransfType,
									LanczosSolverType,VectorWithOffsetType> TimeVectorsBaseType;
			typedef TimeVectorsKrylov<TargettingParamsType,ModelType,WaveFunctionTransfType,
									  LanczosSolverType,VectorWithOffsetType> TimeVectorsKrylovType;
			typedef TimeVectorsRungeKutta<TargettingParamsType,ModelType,WaveFunctionTransfType,
										  LanczosSolverType,VectorWithOffsetType> TimeVectorsRungeKuttaType;
			typedef TimeVectorsSuzukiTrotter<TargettingParamsType,ModelType,WaveFunctionTransfType,
											 LanczosSolverType,VectorWithOffsetType> TimeVectorsSuzukiTrotterType;

			enum {DISABLED,WFT_NOADVANCE,WFT_ADVANCE,COLLAPSE};
			enum {EXPAND_ENVIRON=WaveFunctionTransfType::EXPAND_ENVIRON,
			EXPAND_SYSTEM=WaveFunctionTransfType::EXPAND_SYSTEM,
			INFINITE=WaveFunctionTransfType::INFINITE};
			const static SizeType SYSTEM = ProgramGlobals::SYSTEM;

			static const SizeType parallelRank_ = 0; // Metts needs to support concurrency FIXME

			MettsTargetting(const LeftRightSuperType& lrs,
	 		                const ModelType& model,
			                const TargettingParamsType& mettsStruct,
			                const WaveFunctionTransfType& wft,
			                const SizeType& quantumSector)
			: stage_(DISABLED),
			  lrs_(lrs),
			  model_(model),
			  mettsStruct_(mettsStruct),
			  wft_(wft),
			  quantumSector_(quantumSector),
			  progress_("MettsTargetting",0),
			  currentBeta_(0),
			  applyOpLocal_(lrs),
			  mettsStochastics_(model,mettsStruct.rngSeed),
			  mettsCollapse_(mettsStochastics_,lrs_,mettsStruct),
			  timesWithoutAdvancement_(0),
			  timeVectorsBase_(0),
			  prevDirection_(INFINITE),
			  systemPrev_(),
			  environPrev_()
			{
				if (!wft.isEnabled()) throw PsimagLite::RuntimeError(" MettsTargetting "
							"needs an enabled wft\n");

				RealType tau =mettsStruct_.tau/(mettsStruct_.timeSteps-1);
				SizeType n1 = mettsStruct_.timeSteps;
				SizeType n = mettsStruct_.timeSteps + 1;

				betas_.resize(n1);
				weight_.resize(n);
				targetVectors_.resize(n);

				gsWeight_= 0.0;
				RealType factor = (1.0 - gsWeight_)/(n1+4);
				RealType sum = setOneInterval(factor,PairType(0,n1),tau);
				weight_[n-1] = 2*factor;
				sum += weight_[n-1];

				sum += gsWeight_;
				assert(fabs(sum-1.0)<1e-5);

				PsimagLite::String s (__FILE__);
				s += " Unknown algorithm\n";
				switch (mettsStruct_.algorithm) {
				case TargettingParamsType::KRYLOV:
					timeVectorsBase_ = new TimeVectorsKrylovType(
								currentBeta_,mettsStruct_,betas_,targetVectors_,model_,wft_,lrs_,0);
					break;
				case TargettingParamsType::RUNGE_KUTTA:
					timeVectorsBase_ = new TimeVectorsRungeKuttaType(
								currentBeta_,mettsStruct_,betas_,targetVectors_,model_,wft_,lrs_,0);
					break;
				case TargettingParamsType::SUZUKI_TROTTER:
					timeVectorsBase_ = new TimeVectorsSuzukiTrotterType(
								currentBeta_,mettsStruct_,betas_,targetVectors_,model_,wft_,lrs_,0,0);
					break;
				default:
					throw PsimagLite::RuntimeError(s.c_str());
				}
			}

			~MettsTargetting()
			{
				if (timeVectorsBase_)
					delete timeVectorsBase_;
			}

			const ModelType& model() const { return model_; }

			RealType weight(SizeType i) const
			{
				return weight_[i]; //(allStages(DISABLED)) ? 0.5 : weight_[i];
			}

			RealType gsWeight() const
			{
				return gsWeight_; //(allStages(DISABLED)) ?  0.5 : gsWeight_;
			}

			RealType normSquared(SizeType i) const
			{
				// call to mult will conjugate one of the vectors
				return multiply(targetVectors_[i],targetVectors_[i]); 
			}

			template<typename SomeBasisType>
			void setGs(const typename PsimagLite::Vector<TargetVectorType>::Type& v,
				   const SomeBasisType& someBasis)
			{
				//psi_.set(v,someBasis);
			}

			const RealType& operator[](SizeType i) const
			{
				PsimagLite::String s("MettsTargetting: invalid const operator[]\n");
				throw PsimagLite::RuntimeError(s.c_str());
			}

			RealType& operator[](SizeType i)
			{
				PsimagLite::String s("MettsTargetting: invalid operator[]\n");
				throw PsimagLite::RuntimeError(s.c_str());
			}

			const VectorWithOffsetType& gs() const 
			{
				return targetVectors_[0];
			}

			bool includeGroundStage() const {return (fabs(gsWeight_)>1e-6); }

			SizeType size() const
			{
				if (allStages(DISABLED)) return 1;
				SizeType n = targetVectors_.size();
				if (targetVectors_[n-1].size()==0) n--;
				return n;
			}

			const VectorWithOffsetType& operator()(SizeType i) const
			{
				return targetVectors_[i];
			}

			const LeftRightSuperType& leftRightSuper() const
			{
				return lrs_;
			}

			void evolve(RealType Eg,
			            SizeType direction,
			            const BlockType& block1,
			            const BlockType& block2,
			            SizeType loopNumber)
			{
				typename PsimagLite::Vector<SizeType>::Type sites;
				if (direction==INFINITE)
					utils::blockUnion(sites,block1,block2);
				else sites = block1;

				SizeType n1 = mettsStruct_.timeSteps;

				if (direction==INFINITE) {
					updateStochastics(block1,block2);
					getNewPures(block1,block2);
					return;
				}

				SizeType max = (noStageIs(DISABLED)) ? 1 : n1;
				// Advance or wft each target vector for beta/2
				for (SizeType i=0;i<max;i++) {
					evolve(i,0,n1-1,Eg,direction,sites,loopNumber);
				}

				// compute imag. time evolution:
				calcTimeVectors(PairType(0,n1),Eg,direction);

				// Advance or wft  collapsed vector
				if (targetVectors_[n1].size()>0)
					evolve(n1,n1,n1-1,Eg,direction,sites,loopNumber);

				for (SizeType i=0;i<targetVectors_.size();i++)
					assert(targetVectors_[i].size()==0 || targetVectors_[i].size()==lrs_.super().permutationVector().size());

				cocoon(direction,sites,false);

				if (direction!=INFINITE) {
					if (isAtBorder(direction,sites[0]))
						cocoon(direction,sites,true);
				}

				printEnergies(); // in-situ

				if (stage_!=COLLAPSE) return;

				// collapse
				bool hasCollapsed = mettsCollapse_(targetVectors_[n1],targetVectors_[n1-1],sites,direction);

				if (hasCollapsed) {
					PsimagLite::String s = "  COLLAPSEHERE  ";
					for (SizeType i=0;i<sites.size();i++) {
						test(targetVectors_[n1],targetVectors_[n1],direction,s,sites[i],false);
						if (isAtBorder(direction,sites[0]))
							test(targetVectors_[n1],targetVectors_[n1],direction,s,sites[i],true);
					}
				}
			}

			void load(const PsimagLite::String& f)
			{
				stage_ = WFT_NOADVANCE;
 
 				typename IoType::In io(f);
 
 				TimeSerializerType ts(io,IoType::In::LAST_INSTANCE);
 				for (SizeType i=0;i<targetVectors_.size();i++) targetVectors_[i] = ts.vector(i);
 				currentBeta_ = ts.time();
			}

			void print(std::ostream& os) const
			{
				os<<"MettsWeightsTimeVectors=";
				for (SizeType i=0;i<weight_.size();i++)
					os<<weight_[i]<<" ";
				os<<"\n";
				os<<"MettsWeightGroundState="<<gsWeight_<<"\n";
			}

			void initialGuess(VectorWithOffsetType& v,
			                  const typename PsimagLite::Vector<SizeType>::Type& block) const
			{
				PsimagLite::String s("MettsTargetting: Invalid call to initialGuess\n");
				throw PsimagLite::RuntimeError(s.c_str());
			}

			template<typename IoOutputType>
			void save(const typename PsimagLite::Vector<SizeType>::Type& block,IoOutputType& io) const
			{
				PsimagLite::OstringStream msg;
 				msg<<"Saving state...";
 				progress_.printline(msg,std::cout);

				SizeType marker = 0;
				if (noStageIs(DISABLED)) marker = 1;
				typename PsimagLite::Vector<VectorWithOffsetType>::Type targetVectors(targetVectors_.size());
				if (mettsStruct_.beta>currentBeta_) {
					for (SizeType i=0;i<targetVectors.size();i++)
						targetVectors[i].resize(0);
				} else {
					targetVectors = targetVectors_;
				}
				TimeSerializerType ts(currentBeta_,block[0],targetVectors,marker);
 				ts.save(io);
			}

			RealType time() const { return 0; }

			void updateOnSiteForTimeDep(BasisWithOperatorsType& basisWithOps) const
			{
				// nothing to do here
			}

		private:

			void evolve(SizeType index,
			            SizeType start,
			            SizeType indexAdvance,
			            RealType Eg,
			            SizeType direction,
			            const typename PsimagLite::Vector<SizeType>::Type& block,
			            SizeType loopNumber)
			{
				if (index==0 && start==0)
					advanceCounterAndComputeStage(block);

				PsimagLite::OstringStream msg;
				msg<<"Evolving, stage="<<getStage()<<" loopNumber="<<loopNumber;
				msg<<" Eg="<<Eg;
				progress_.printline(msg,std::cout);
				advanceOrWft(index,indexAdvance,direction,block);
			}

			void advanceCounterAndComputeStage(const typename PsimagLite::Vector<SizeType>::Type& block)
			{
				if (stage_!=COLLAPSE) stage_=WFT_NOADVANCE;

				if (stage_==COLLAPSE) {
					if (!allSitesCollapsed()) {
						if (sitesCollapsed_.size()>2*model_.geometry().numberOfSites())
							throw PsimagLite::RuntimeError("advanceCounterAndComputeStage\n");
						printAdvancement();
						return;
					}

					sitesCollapsed_.clear();
					stage_ = WFT_NOADVANCE;
					timesWithoutAdvancement_=0;
					currentBeta_ = 0;
					PsimagLite::OstringStream msg;
					SizeType n1 = mettsStruct_.timeSteps;
					RealType x = std::norm(targetVectors_[n1]);
					msg<<"Changing direction, setting collapsed with norm="<<x;
					progress_.printline(msg,std::cout);
					for (SizeType i=0;i<n1;i++)
						targetVectors_[i] = targetVectors_[n1];
					timeVectorsBase_->timeHasAdvanced();
					printAdvancement();
					return;
				}

				if (timesWithoutAdvancement_ < mettsStruct_.advanceEach) {
					timesWithoutAdvancement_++;
					printAdvancement();
					return;
				}

				if (stage_!=COLLAPSE && currentBeta_<mettsStruct_.beta) {
					stage_ = WFT_ADVANCE;
					currentBeta_ += mettsStruct_.tau;
					timesWithoutAdvancement_=0;
					printAdvancement();
					return;
				}

				if (stage_!=COLLAPSE && currentBeta_>=mettsStruct_.beta && block[0]!=block.size()) {
					printAdvancement();
					return;
				}

				if (stage_!=COLLAPSE && currentBeta_>=mettsStruct_.beta) {
					stage_ = COLLAPSE;
					sitesCollapsed_.clear();
					SizeType n1 = mettsStruct_.timeSteps;
					targetVectors_[n1].resize(0);
					timesWithoutAdvancement_=0;
					printAdvancement();
					return;
				}
			}

			void printAdvancement() const
			{
				PsimagLite::OstringStream msg2;
				msg2<<"Steps without advance: "<<timesWithoutAdvancement_;
				if (timesWithoutAdvancement_>0)
					progress_.printline(msg2,std::cout);
			}

			void advanceOrWft(SizeType index,
			                  SizeType indexAdvance,
			                  SizeType systemOrEnviron,
			                  const typename PsimagLite::Vector<SizeType>::Type& block)
			{
				if (targetVectors_[index].size()==0) return;
				assert(std::norm(targetVectors_[index])>1e-6);
				typename PsimagLite::Vector<SizeType>::Type nk;
				mettsCollapse_.setNk(nk,block);

				if (stage_== WFT_NOADVANCE || stage_== WFT_ADVANCE || stage_==COLLAPSE) {
					SizeType advance = index;
					if (stage_ == WFT_ADVANCE) {
						advance = indexAdvance;
						timeVectorsBase_->timeHasAdvanced();
					}
					// don't advance the collapsed vector because we'll recompute
					if (index==weight_.size()-1) advance=index;
					PsimagLite::OstringStream msg;
					msg<<"I'm calling the WFT now";
					progress_.printline(msg,std::cout);

					VectorWithOffsetType phiNew; // same sectors as g.s.
					//phiNew.populateSectors(lrs_.super());
					assert(std::norm(targetVectors_[advance])>1e-6);

					phiNew.populateSectors(lrs_.super());
					// OK, now that we got the partition number right, let's wft:
					wft_.setInitialVector(phiNew,targetVectors_[advance],lrs_,nk);
					phiNew.collapseSectors();
					assert(std::norm(phiNew)>1e-6);
					targetVectors_[index] = phiNew;
				} else {
					assert(false);
				}
			}

			void updateStochastics(const typename PsimagLite::Vector<SizeType>::Type& block1,const typename PsimagLite::Vector<SizeType>::Type& block2)
			{
				SizeType linSize = model_.geometry().numberOfSites();
				typename PsimagLite::Vector<SizeType>::Type tqn(2,0);
				if (model_.params().targetQuantumNumbers.size()>=2) {
					tqn[0] = SizeType(round(model_.params().targetQuantumNumbers[0]*linSize));
					tqn[1] = SizeType(round(model_.params().targetQuantumNumbers[1]*linSize));
				} else {
					tqn[0] = model_.params().electronsUp;
					tqn[1] = model_.params().electronsDown;
				}
				SizeType qn = BasisType::pseudoQuantumNumber(tqn);
				mettsStochastics_.update(qn,block1,block2,mettsStruct_.rngSeed);
			}

			SizeType getPartition() const
			{
				SizeType total = lrs_.super().partition()-1;
				for (SizeType i=0;i<total;i++) {
					// Do only one sector unless doing su(2) with j>0, then do all m's
					if (lrs_.super().pseudoEffectiveNumber(
						lrs_.super().partition(i))==quantumSector_ )
						return i;
				}
				throw PsimagLite::RuntimeError("MettsTargetting: getPartition()\n");
			}

			// direction here is INFINITE
			void getNewPures(const typename PsimagLite::Vector<SizeType>::Type& block1,const typename PsimagLite::Vector<SizeType>::Type& block2)
			{
				typename PsimagLite::Vector<SizeType>::Type alphaFixed(block1.size());
				for (SizeType i=0;i<alphaFixed.size();i++)
					alphaFixed[i] = mettsStochastics_.chooseRandomState(block1[i]);

				typename PsimagLite::Vector<SizeType>::Type betaFixed(block2.size());
				for (SizeType i=0;i<betaFixed.size();i++)
					betaFixed[i] = mettsStochastics_.chooseRandomState(block2[i]);

				PsimagLite::OstringStream msg;
				msg<<"New pures for ";
				for (SizeType i=0;i<alphaFixed.size();i++)
					msg<<" site="<<block1[i]<<" is "<<alphaFixed[i];
				msg<<" and for ";
				for (SizeType i=0;i<betaFixed.size();i++)
					msg<<" site="<<block2[i]<<" is "<<betaFixed[i];
				progress_.printline(msg,std::cerr);

				const MatrixType& transformSystem =  wft_.transform(ProgramGlobals::SYSTEM);
				VectorType newVector1(transformSystem.n_row(),0);

				typename PsimagLite::Vector<SizeType>::Type nk1;
				mettsCollapse_.setNk(nk1,block1);
				SizeType alphaFixedVolume = mettsCollapse_.volumeOf(alphaFixed,nk1);

				getNewPure(newVector1,pureVectors_.first,ProgramGlobals::SYSTEM,
				           alphaFixedVolume,lrs_.left(),transformSystem,block1);
				pureVectors_.first = newVector1;

				const MatrixType& transformEnviron = 
				                        wft_.transform(ProgramGlobals::ENVIRON);
				VectorType newVector2(transformEnviron.n_row(),0);

				typename PsimagLite::Vector<SizeType>::Type nk2;
				mettsCollapse_.setNk(nk2,block2);
				SizeType betaFixedVolume = mettsCollapse_.volumeOf(betaFixed,nk2);
				getNewPure(newVector2,pureVectors_.second,ProgramGlobals::ENVIRON,
						   betaFixedVolume,lrs_.right(),transformEnviron,block2);
				pureVectors_.second = newVector2;
				setFromInfinite(targetVectors_[0],lrs_);
				assert(std::norm(targetVectors_[0])>1e-6);

				systemPrev_.fixed = alphaFixedVolume;
				systemPrev_.permutationInverse = lrs_.left().permutationInverse();
				environPrev_.fixed = betaFixedVolume;
				environPrev_.permutationInverse = lrs_.right().permutationInverse();
			}

			void getFullVector(typename PsimagLite::Vector<RealType>::Type& v,SizeType m,const LeftRightSuperType& lrs) const
			{
				int offset = lrs.super().partition(m);
				int total = lrs.super().partition(m+1) - offset;

				PackIndicesType pack(lrs.left().size());
				v.resize(total);
				assert(PsimagLite::norm(pureVectors_.first)>1e-6);
				assert(PsimagLite::norm(pureVectors_.second)>1e-6);
				for (int i=0;i<total;i++) {
					SizeType alpha,beta;
					pack.unpack(alpha,beta,lrs.super().permutation(i+offset));
					v[i] = pureVectors_.first[alpha] * pureVectors_.second[beta];
				}
			}

			void getNewPure(VectorType& newVector,
			                VectorType& oldVector,
			                SizeType direction,
			                SizeType alphaFixed,
			                const BasisWithOperatorsType& basis,
			                const MatrixType& transform,
			                const typename PsimagLite::Vector<SizeType>::Type& block)
			{
				if (oldVector.size()==0)
					setInitialPure(oldVector,block);
				VectorType tmpVector;
				if (transform.n_row()==0) {
					tmpVector = oldVector;
					assert(PsimagLite::norm(tmpVector)>1e-6);
				} else {
					delayedTransform(tmpVector,oldVector,direction,transform,block);
					assert(PsimagLite::norm(tmpVector)>1e-6);
				}
				SizeType ns = tmpVector.size();
				typename PsimagLite::Vector<SizeType>::Type nk;
				mettsCollapse_.setNk(nk,block);
				SizeType volumeOfNk = mettsCollapse_.volumeOf(nk);
				SizeType newSize =  (transform.n_col()==0) ? (ns*ns) : 
							transform.n_col() * volumeOfNk;
				newVector.resize(newSize);
				for (SizeType alpha=0;alpha<newVector.size();alpha++)
					newVector[alpha] = 0;

				for (SizeType alpha=0;alpha<ns;alpha++) {
					SizeType gamma = (direction==ProgramGlobals::SYSTEM) ?
						basis.permutationInverse(alpha + alphaFixed*ns) :
						basis.permutationInverse(alphaFixed + alpha*volumeOfNk);
					newVector[gamma] = tmpVector[alpha];
				}
//				for (SizeType gamma=0;gamma<newVector.size();gamma++) {
//					newVector[gamma] = 0;
//					for (SizeType alpha=0;alpha<ns;alpha++) {
//						SizeType gammaPrime = (direction==ProgramGlobals::SYSTEM) ?
//									basis.permutationInverse(alpha + alphaFixed*ns) :
//									basis.permutationInverse(alphaFixed + alpha*nk);

//						if (gamma == gammaPrime)
//							newVector[gamma] += tmpVector[alpha];
//					}
//				}
				PsimagLite::OstringStream msg2;
				msg2<<"Old size of pure is "<<ns<<" norm="<<PsimagLite::norm(tmpVector);
				progress_.printline(msg2,std::cerr);
				PsimagLite::OstringStream msg;
				msg<<"New size of pure is "<<newSize<<" norm="<<PsimagLite::norm(newVector);
				progress_.printline(msg,std::cerr);
				assert(PsimagLite::norm(newVector)>1e-6);
			}

			void delayedTransform(VectorType& newVector,
								  VectorType& oldVector,
								  SizeType direction,
								  const MatrixType& transform,
								  const typename PsimagLite::Vector<SizeType>::Type& block)
			{
				assert(oldVector.size()==transform.n_row());

				typename PsimagLite::Vector<SizeType>::Type nk;
				mettsCollapse_.setNk(nk,block);
				SizeType ne = mettsCollapse_.volumeOf(nk);
				
				const typename PsimagLite::Vector<SizeType>::Type& permutationInverse = (direction==SYSTEM)
				? systemPrev_.permutationInverse : environPrev_.permutationInverse;
				SizeType nsPrev = permutationInverse.size()/ne;
				
				newVector.resize(transform.n_col());
				//newVector = oldVector * transform;
				for (SizeType gamma=0;gamma<newVector.size();gamma++) {
					newVector[gamma] = 0;
					for (SizeType alpha=0;alpha<nsPrev;alpha++) {
						SizeType noPermIndex =  (direction==SYSTEM)
								   ? alpha + systemPrev_.fixed*nsPrev
								   : environPrev_.fixed + alpha*ne;
						
						SizeType gammaPrime = permutationInverse[noPermIndex];
						
						assert(gammaPrime<transform.n_row());
						newVector[gamma] += transform(gammaPrime,gamma) *
									  oldVector[gammaPrime];
					}
				}
			}

			void simpleTransform(VectorType& newVector,
								  VectorType& oldVector,
								  const MatrixType& transform)
			{
				assert(oldVector.size()==transform.n_row());


				for (SizeType gamma=0;gamma<newVector.size();gamma++) {
					newVector[gamma] = 0;
					for (SizeType gammaPrime=0;gammaPrime<oldVector.size();gammaPrime++) {
						newVector[gamma] += transform(gammaPrime,gamma) *
									  oldVector[gammaPrime];
					}
				}
			}

			void setInitialPure(VectorType& oldVector,const typename PsimagLite::Vector<SizeType>::Type& block)
			{
				int offset = (block[0]==block.size()) ? -block.size() : block.size();
				typename PsimagLite::Vector<SizeType>::Type blockCorrected = block;
				for (SizeType i=0;i<blockCorrected.size();i++)
					blockCorrected[i] += offset;

				typename PsimagLite::Vector<SizeType>::Type nk;
				mettsCollapse_.setNk(nk,blockCorrected);
				SizeType volumeOfNk = mettsCollapse_.volumeOf(nk);
				typename PsimagLite::Vector<SizeType>::Type alphaFixed(nk.size());
				for (SizeType i=0;i<alphaFixed.size();i++)
					alphaFixed[i] = mettsStochastics_.chooseRandomState(blockCorrected[i]);

				PsimagLite::OstringStream msg;
				msg<<"New pures for site ";
				for (SizeType i=0;i<blockCorrected.size();i++)
					msg<<blockCorrected[i]<<" ";
				msg<<" is "<<alphaFixed;
				progress_.printline(msg,std::cerr);

				SizeType volumeOfAlphaFixed = mettsCollapse_.volumeOf(alphaFixed,nk);

				oldVector.resize(volumeOfNk);
				assert(volumeOfAlphaFixed<oldVector.size());
				for (SizeType i=0;i<oldVector.size();i++) {
					oldVector[i] = (i==volumeOfAlphaFixed) ? 1 : 0;
				}
				assert(PsimagLite::norm(oldVector)>1e-6);
			}

			void setFromInfinite(VectorWithOffsetType& phi,const LeftRightSuperType& lrs) const
			{
				phi.populateSectors(lrs.super());
				for (SizeType ii=0;ii<phi.sectors();ii++) {
					SizeType i0 = phi.sector(ii);
					VectorType v;
					getFullVector(v,i0,lrs);
					RealType tmpNorm = PsimagLite::norm(v);
					if (fabs(tmpNorm-1.0)<1e-6) {
						SizeType j = lrs.super().qn(lrs.super().partition(i0));
						typename PsimagLite::Vector<SizeType>::Type qns = BasisType::decodeQuantumNumber(j);
						std::cerr<<"setFromInfinite: qns= ";
						for (SizeType k=0;k<qns.size();k++) std::cerr<<qns[k]<<" ";
						std::cerr<<"\n";
					}
					phi.setDataInSector(v,i0);
				}
				phi.collapseSectors();
				assert(std::norm(phi)>1e-6);
			}

			// in situ computation:
			void cocoon(SizeType direction,const typename PsimagLite::Vector<SizeType>::Type& block,bool corner)
			{
				for (SizeType i=0;i<block.size();i++)
					cocoon(direction,block[i],corner);
			}

			// in situ computation:
			void cocoon(SizeType direction,SizeType site,bool corner)
			{
				std::cerr<<"-------------&*&*&* In-situ measurements start\n";
				
				for (SizeType j=0;j<targetVectors_.size();j++) {
					PsimagLite::String s = "<P"+ttos(j)+"|A|P"+ttos(j)+">";
					SizeType site2 = test(targetVectors_[j],targetVectors_[j],direction,s,site,corner);
					if (stage_==COLLAPSE && j==0) sitesCollapsed_.push_back(site2);
				}
				std::cerr<<"-------------&*&*&* In-situ measurements end\n";
			}

			void checkOrder(SizeType i) const
			{
				if (i==0) return;
				for (SizeType j=0;j<i;j++) {
					if (stage_ == DISABLED) {
						PsimagLite::String s ="TST:: Seeing tst site "+ttos(mettsStruct_.sites[i]);
						s =s + " before having seen";
						s = s + " site "+ttos(j);
						s = s +". Please order your tst sites in order of appearance.\n";
						throw PsimagLite::RuntimeError(s);
					}
				}
			}

			bool allStages(SizeType x) const
			{
				if (stage_!=x) return false;
				return true;
			}

			bool noStageIs(SizeType x) const
			{
				//for (SizeType i=0;i<stage_.size();i++)
				if (stage_==x) return false;
				return true;
			}

			PsimagLite::String getStage() const
			{
				switch (stage_) {
					case DISABLED:
						return "Disabled";
						break;
					case COLLAPSE:
						return "Collapsing";
						break;
					case WFT_ADVANCE:
						return "WFT with time stepping";
						break;
					case WFT_NOADVANCE:
						return "WFT without time change";
						break;
				}
				return "undefined";
			}

			void calcTimeVectors(const PairType& startEnd,
			                     RealType Eg,
			                     SizeType systemOrEnviron)
			{
				const VectorWithOffsetType& phi = targetVectors_[startEnd.first];
				PsimagLite::OstringStream msg;
				msg<<" vector number "<<startEnd.first<<" has norm ";
				msg<<std::norm(phi);
				progress_.printline(msg,std::cout);
				if (std::norm(phi)<1e-6) setFromInfinite(targetVectors_[startEnd.first],lrs_);
				bool allOperatorsApplied = (noStageIs(DISABLED));
				timeVectorsBase_->calcTimeVectors(startEnd,Eg,phi,systemOrEnviron,allOperatorsApplied);
				normalizeVectors(startEnd);
			}

			void normalizeVectors(const PairType& startEnd)
			{
				for (SizeType i=startEnd.first+1;i<startEnd.second;i++) {
					VectorWithOffsetType v = targetVectors_[i];
					RealType x = 1.0/std::norm(v);
					targetVectors_[i]= x* v;
				}
			}

			RealType setOneInterval(const RealType& factor,
			                        const PairType& startEnd,
			                        const RealType& tau)
			{
				RealType sum = 0;
				for (SizeType i=startEnd.first;i<startEnd.second;i++) {
					betas_[i] = (i-startEnd.first)*tau;
					weight_[i] = factor;
					sum += weight_[i];
				}
				
				sum -= weight_[startEnd.first];
				sum -= weight_[startEnd.second-1];
				weight_[startEnd.first] = weight_[startEnd.second-1] = 2*factor;
				sum += weight_[startEnd.second-1];
				sum += weight_[startEnd.first];
				return sum;
			}

			void findElectronsOfOneSite(typename PsimagLite::Vector<SizeType>::Type& electrons,SizeType site) const
			{
				typename PsimagLite::Vector<SizeType>::Type block(1,site);
				typename ModelType::HilbertBasisType basis;
				typename PsimagLite::Vector<SizeType>::Type quantumNumbs;
				model_.setNaturalBasis(basis,quantumNumbs,block);
				model_.findElectrons(electrons,basis,site);
			}

			SizeType test(const VectorWithOffsetType& src1,
			            const VectorWithOffsetType& src2,
			            SizeType systemOrEnviron,
			            const PsimagLite::String& label,
			            SizeType site,
			            bool corner) const
			{
				VectorWithOffsetType dest;
				OperatorType A = getObservableToTest(model_.params().model);
				typename PsimagLite::Vector<SizeType>::Type electrons;
				findElectronsOfOneSite(electrons,site);
				FermionSign fs(lrs_.left(),electrons);
				applyOpLocal_(dest,src1,A,fs,systemOrEnviron,corner);

				RealType sum = 0;
				for (SizeType ii=0;ii<dest.sectors();ii++) {
					SizeType i = dest.sector(ii);
					SizeType offset1 = dest.offset(i);
					for (SizeType jj=0;jj<src2.sectors();jj++) {
						SizeType j = src2.sector(jj);
						SizeType offset2 = src2.offset(j);
						if (i!=j) continue; //throw PsimagLite::RuntimeError("Not same sector\n");
						for (SizeType k=0;k<dest.effectiveSize(i);k++) 
							sum+= dest[k+offset1] * std::conj(src2[k+offset2]);
					}
				}
				RealType nor = std::norm(src1);
				SizeType site2 = site;
				if (corner) {
					if (site2==1) {
						site2=0;
					} else {
						site2 = model_.geometry().numberOfSites() -1;
					}
				}

				std::cerr<<site2<<" "<<sum<<" "<<" "<<currentBeta_;
				std::cerr<<" "<<label<<" "<<nor<<" "<<std::norm(src2);
				std::cerr<<" "<<std::norm(dest)<<"    "<<sum/(nor*nor)<<"\n";

				return site2;
			}

			OperatorType getObservableToTest(const PsimagLite::String& modelName) const
			{
				OperatorType A;
				SizeType site = 0; // sites.first; <-- site-dependent Hilbert space not supported by METTS

				if (modelName=="HubbardOneBand") {
					PsimagLite::CrsMatrix<RealType> tmpC(model_.naturalOperator("nup",site,0));
					A.data = tmpC;
					A.fermionSign = 1;
					return A;
				}
				if (modelName=="FeAsBasedSc" || modelName=="FeAsBasedScExtended") {
					PsimagLite::CrsMatrix<RealType> tmpC(model_.naturalOperator("c",site,0));
					PsimagLite::CrsMatrix<RealType> tmpCdagger;
					transposeConjugate(tmpCdagger,tmpC);
					multiply(A.data,tmpCdagger,tmpC);
					A.fermionSign = 1;
					return A;
				}
				PsimagLite::String s(__FILE__);
				s += " " + ttos(__LINE__) + "\n";
				s += "Model " + modelName + " not supported by MettsTargetting\n";
				throw PsimagLite::RuntimeError(s.c_str());
			}

			void printEnergies() const
			{
				for (SizeType i=0;i<targetVectors_.size();i++)
					printEnergies(targetVectors_[i],i);
			}

			void printEnergies(const VectorWithOffsetType& phi,SizeType whatTarget) const
			{
				for (SizeType ii=0;ii<phi.sectors();ii++) {
					SizeType i = phi.sector(ii);
					printEnergies(phi,whatTarget,i);
				}
			}

			void printEnergies(const VectorWithOffsetType& phi,SizeType whatTarget, SizeType i0) const
			{
				SizeType p = lrs_.super().findPartitionNumber(phi.offset(i0));
				typename ModelType::ModelHelperType modelHelper(p,lrs_);
						//,useReflection_);
				typename LanczosSolverType::LanczosMatrixType lanczosHelper(&model_,&modelHelper);


				SizeType total = phi.effectiveSize(i0);
				TargetVectorType phi2(total);
				phi.extract(phi2,i0);
				TargetVectorType x(total);
				lanczosHelper.matrixVectorProduct(x,phi2);
				PsimagLite::OstringStream msg;
				msg<<"Hamiltonian average at beta="<<currentBeta_<<" for target="<<whatTarget;
				msg<<" sector="<<i0<<" <phi(t)|H|phi(t)>="<<(phi2*x)<<" <phi(t)|phi(t)>="<<(phi2*phi2);
				progress_.printline(msg,std::cout);
			}

			bool isAtBorder(SizeType direction,SizeType site) const
			{
				if (direction==EXPAND_SYSTEM && site+2==model_.geometry().numberOfSites())
					return true;
				if (direction!=EXPAND_SYSTEM && site==1)
					return true;
				return false;
			}

			bool allSitesCollapsed() const
			{
				SizeType n = model_.geometry().numberOfSites();
				for (SizeType i=0;i<n;i++) {
					bool seen = (std::find(sitesCollapsed_.begin(),sitesCollapsed_.end(),i) != sitesCollapsed_.end());
					if (!seen) return false;
				}
				return true;
			}

			SizeType stage_;
			const LeftRightSuperType& lrs_;
			const ModelType& model_;
			const TargettingParamsType& mettsStruct_;
			const WaveFunctionTransfType& wft_;
			const SizeType& quantumSector_;
			PsimagLite::ProgressIndicator progress_;
			RealType currentBeta_;
			typename PsimagLite::Vector<RealType>::Type betas_,weight_;
			typename PsimagLite::Vector<VectorWithOffsetType>::Type targetVectors_;
			RealType gsWeight_;
			ApplyOperatorType applyOpLocal_;
			MettsStochasticsType mettsStochastics_;
			MettsCollapseType mettsCollapse_;
			SizeType timesWithoutAdvancement_;
			TimeVectorsBaseType* timeVectorsBase_;
			SizeType prevDirection_;
			MettsPrev systemPrev_;
			MettsPrev environPrev_;
			std::pair<VectorType,VectorType> pureVectors_;
			typename PsimagLite::Vector<SizeType>::Type sitesCollapsed_;
	};     //class MettsTargetting

	template<
		template<typename,typename,typename> class LanczosSolverTemplate,
			template<typename,typename> class InternalProductTemplate,
 		template<typename,typename> class WaveFunctionTransfTemplate,
			typename ModelType_,
			typename IoType_,
   			template<typename> class VectorWithOffsetTemplate>
	std::ostream& operator<<(std::ostream& os,
			const MettsTargetting<LanczosSolverTemplate,
			InternalProductTemplate,
			WaveFunctionTransfTemplate,ModelType_,IoType_,
			VectorWithOffsetTemplate>& tst)
	{
		tst.print(os);
		return os;
	}
} // namespace Dmrg

#endif //METTS_TARGETTING_H
