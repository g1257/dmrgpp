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

namespace Dmrg {
	template<
			template<typename,typename,typename> class LanczosSolverTemplate,
			template<typename,typename> class InternalProductTemplate,
	 		template<typename,typename> class WaveFunctionTransfTemplate,
    			typename ModelType_,
    			typename ConcurrencyType_,
    			typename IoType_,
       			template<typename> class VectorWithOffsetTemplate>
	class MettsTargetting  {

			struct MettsPrev {
				MettsPrev() : fixed(0),permutationInverse(0) { }
				size_t fixed;
				std::vector<size_t> permutationInverse;
			};

		public:
			typedef ModelType_ ModelType;
			typedef ConcurrencyType_ ConcurrencyType;
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
			typedef std::vector<RealType> VectorType;
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

			enum {DISABLED,WFT_NOADVANCE,WFT_ADVANCE,COLLAPSE};
			enum {EXPAND_ENVIRON=WaveFunctionTransfType::EXPAND_ENVIRON,
			EXPAND_SYSTEM=WaveFunctionTransfType::EXPAND_SYSTEM,
			INFINITE=WaveFunctionTransfType::INFINITE};
			const static size_t SYSTEM = ProgramGlobals::SYSTEM;

			static const size_t parallelRank_ = 0; // Metts needs to support concurrency FIXME

			MettsTargetting(const LeftRightSuperType& lrs,
	 		                const ModelType& model,
			                const TargettingParamsType& mettsStruct,
			                const WaveFunctionTransfType& wft,
			                const size_t& quantumSector)
			: stage_(DISABLED),
			  lrs_(lrs),
			  model_(model),
			  mettsStruct_(mettsStruct),
			  wft_(wft),
			  quantumSector_(quantumSector),
			  progress_("MettsTargetting(AlphaStage)",0),
			  currentBeta_(0),
			  applyOpLocal_(lrs),
			  mettsStochastics_(model,mettsStruct.rngSeed),
			  mettsCollapse_(mettsStochastics_,lrs_,mettsStruct),
			  timesWithoutAdvancement_(0),
			  prevDirection_(INFINITE),
			  systemPrev_(),
			  environPrev_()
			{
				if (!wft.isEnabled()) throw std::runtime_error(" MettsTargetting "
							"needs an enabled wft\n");

				RealType tau =mettsStruct_.tau/(mettsStruct_.timeSteps-1);
				size_t n1 = mettsStruct_.timeSteps;
				size_t n = mettsStruct_.timeSteps + 1;

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
			}

			const ModelType& model() const { return model_; }

			RealType weight(size_t i) const
			{
				return weight_[i]; //(allStages(DISABLED)) ? 0.5 : weight_[i];
			}

			RealType gsWeight() const
			{
				return gsWeight_; //(allStages(DISABLED)) ?  0.5 : gsWeight_;
			}

			RealType normSquared(size_t i) const
			{
				// call to mult will conjugate one of the vectors
				return multiply(targetVectors_[i],targetVectors_[i]); 
			}

			template<typename SomeBasisType>
			void setGs(const std::vector<TargetVectorType>& v,
				   const SomeBasisType& someBasis)
			{
				//psi_.set(v,someBasis);
			}

			const RealType& operator[](size_t i) const
			{
				std::string s("MettsTargetting: invalid const operator[]\n");
				throw std::runtime_error(s.c_str());
			}

			RealType& operator[](size_t i)
			{
				std::string s("MettsTargetting: invalid operator[]\n");
				throw std::runtime_error(s.c_str());
			}

			const VectorWithOffsetType& gs() const 
			{
				return targetVectors_[0];
			}

			bool includeGroundStage() const {return (fabs(gsWeight_)>1e-6); }

			size_t size() const
			{
				if (allStages(DISABLED)) return 1;
				size_t n = targetVectors_.size();
				if (targetVectors_[n-1].size()==0) n--;
				return n;
			}

			const VectorWithOffsetType& operator()(size_t i) const
			{
				return targetVectors_[i];
			}

			const LeftRightSuperType& leftRightSuper() const
			{
				return lrs_;
			}

			void evolve(RealType Eg,
			            size_t direction,
			            const BlockType& block1,
			            const BlockType& block2,
			            size_t loopNumber)
			{
				assert(block1.size()==1);

				PairType sites(block1[0],block2[0]);
				size_t n1 = mettsStruct_.timeSteps;
				
				updateStochastics(sites);

				if (direction==INFINITE) {
					getNewPures(sites,n1);
					return;
				}

				// Advance or wft each target vector for beta/2
				for (size_t i=0;i<n1;i++) {
					evolve(i,0,n1-1,Eg,direction,sites,loopNumber);
				}

				// compute imag. time evolution:
				calcTimeVectors(PairType(0,n1),Eg,direction);

				// Advance or wft  collapsed vector
				if (targetVectors_[n1].size()>0)
					evolve(n1,n1,n1-1,Eg,direction,sites,loopNumber);

				for (size_t i=0;i<targetVectors_.size();i++)
					assert(targetVectors_[i].size()==0 || targetVectors_[i].size()==lrs_.super().permutationVector().size());

				cocoon(direction,sites,false);

				if (sites.first==sites.second) {
					if (isAtBorder(direction,sites.first))
						cocoon(direction,sites,true);
				}

				printEnergies(); // in-situ

				if (stage_!=COLLAPSE) return;

				// collapse
				bool hasCollapsed = mettsCollapse_(targetVectors_[n1],targetVectors_[n1-1],sites.first,direction);

				if (hasCollapsed) {
					//throw std::runtime_error("collapsed testing\n");
					std::string s = "  COLLAPSEHERE  ";
					test(targetVectors_[n1],targetVectors_[n1],direction,s,sites,false);
					if (isAtBorder(direction,sites.first))
						test(targetVectors_[n1],targetVectors_[n1],direction,s,sites,true);
				}
			}

			void load(const std::string& f)
			{
				//throw std::runtime_error("Metts: load() unimplemented\n");
				stage_ = WFT_NOADVANCE;
 
 				typename IoType::In io(f);
 
 				TimeSerializerType ts(io,IoType::In::LAST_INSTANCE);
 				for (size_t i=0;i<targetVectors_.size();i++) targetVectors_[i] = ts.vector(i);
 				currentBeta_ = ts.time();
 
				/*int site = 0;
				io.readline(site,"#TCENTRALSITE=",IoType::In::LAST_INSTANCE);
								if (site<0) throw std::runtime_error("Metts::load(...): site cannot be negative\n");
 				targetVectors_[0].load(io,"PSI");
				*/
			}

			void print(std::ostream& os) const
			{
				os<<"MettsWeightsTimeVectors=";
				for (size_t i=0;i<weight_.size();i++)
					os<<weight_[i]<<" ";
				os<<"\n";
				os<<"MettsWeightGroundState="<<gsWeight_<<"\n";
			}

			void initialGuess(VectorWithOffsetType& v,size_t nk) const
			{
				std::string s("MettsTargetting: Invalid call to initialGuess\n");
				throw std::runtime_error(s.c_str());
			}

			template<typename IoOutputType>
			void save(const std::vector<size_t>& block,IoOutputType& io) const
			{
				std::ostringstream msg;
 				msg<<"Saving state...";
 				progress_.printline(msg,std::cout);

				size_t marker = 0;
				if (noStageIs(DISABLED)) marker = 1;
				std::vector<VectorWithOffsetType> targetVectors(targetVectors_.size());
				if (mettsStruct_.beta>currentBeta_) {
					for (size_t i=0;i<targetVectors.size();i++)
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

			void evolve(size_t index,
			            size_t start,
			            size_t indexAdvance,
			            RealType Eg,
			            size_t direction,
			            std::pair<size_t,size_t> sites,
			            size_t loopNumber)
			{
				if (index==0 && start==0) advanceCounterAndComputeStage(sites.first);

				std::ostringstream msg;
				msg<<"Evolving, stage="<<getStage()<<" loopNumber="<<loopNumber;
				msg<<" Eg="<<Eg;
				progress_.printline(msg,std::cout);
				assert(sites.first==sites.second);
				size_t nk = model_.hilbertSize(sites.first);
				advanceOrWft(index,indexAdvance,direction,nk,sites);
			}

			void advanceCounterAndComputeStage(size_t site)
			{
				if (stage_!=COLLAPSE) stage_=WFT_NOADVANCE;

				if (stage_==COLLAPSE) {
					if (!allSitesCollapsed()) {
						if (sitesCollapsed_.size()>2*model_.geometry().numberOfSites())
							throw std::runtime_error("advanceCounterAndComputeStage\n");
						printAdvancement();
						return;
					}

					sitesCollapsed_.clear();
					stage_ = WFT_NOADVANCE;
					timesWithoutAdvancement_=0;
					currentBeta_ = 0;
					std::ostringstream msg;
					size_t n1 = mettsStruct_.timeSteps;
					RealType x = std::norm(targetVectors_[n1]);
					msg<<"Changing direction, setting collapsed with norm="<<x;
					progress_.printline(msg,std::cout);
					targetVectors_[0] = targetVectors_[n1];
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

				if (stage_!=COLLAPSE && currentBeta_>=mettsStruct_.beta && site!=1) {
					printAdvancement();
					return;
				}

				if (stage_!=COLLAPSE && currentBeta_>=mettsStruct_.beta) {
					stage_ = COLLAPSE;
					sitesCollapsed_.clear();
					size_t n1 = mettsStruct_.timeSteps;
					targetVectors_[n1].resize(0);
					timesWithoutAdvancement_=0;
					printAdvancement();
					return;
				}
			}

			void printAdvancement() const
			{
				std::ostringstream msg2;
				msg2<<"Steps without advance: "<<timesWithoutAdvancement_;
				if (timesWithoutAdvancement_>0)
					progress_.printline(msg2,std::cout);
			}

			void advanceOrWft(size_t index,
					  size_t indexAdvance,
			                  size_t systemOrEnviron,
					  size_t nk,
					  std::pair<size_t,size_t> sites)
			{
				if (targetVectors_[index].size()==0) return;
				assert(std::norm(targetVectors_[index])>1e-6);

				if (stage_== WFT_NOADVANCE || stage_== WFT_ADVANCE || stage_==COLLAPSE) {
					size_t advance = index;
					if (stage_ == WFT_ADVANCE) advance = indexAdvance;
					// don't advance the collapsed vector because we'll recompute
					if (index==weight_.size()-1) advance=index;
					std::ostringstream msg;
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

			void updateStochastics(const PairType& sites)
			{
				size_t linSize = model_.geometry().numberOfSites();
				std::vector<size_t> tqn(2,0);
				if (model_.params().targetQuantumNumbers.size()>=2) {
					tqn[0] = size_t(round(model_.params().targetQuantumNumbers[0]*linSize));
					tqn[1] = size_t(round(model_.params().targetQuantumNumbers[1]*linSize));
				} else {
					tqn[0] = model_.params().electronsUp;
					tqn[1] = model_.params().electronsDown;
				}
				size_t qn = BasisType::pseudoQuantumNumber(tqn);
				mettsStochastics_.update(qn,sites,mettsStruct_.rngSeed);
			}

			size_t getPartition() const
			{
				size_t total = lrs_.super().partition()-1;
				for (size_t i=0;i<total;i++) {
					// Do only one sector unless doing su(2) with j>0, then do all m's
					if (lrs_.super().pseudoEffectiveNumber(
						lrs_.super().partition(i))==quantumSector_ )
						return i;
				}
				throw std::runtime_error("MettsTargetting: getPartition()\n");
			}

			void getNewPures(const PairType& sites,size_t n1)
			{
				size_t alphaFixed = mettsStochastics_.chooseRandomState(sites.first);
				size_t betaFixed = mettsStochastics_.chooseRandomState(sites.second);

				std::ostringstream msg;
				msg<<"New pures for site "<<sites.first<<" is "<<alphaFixed;
				msg<<" and for site "<<sites.second<<" is "<<betaFixed;
				progress_.printline(msg,std::cerr);

				const MatrixType& transformSystem =  wft_.transform(ProgramGlobals::SYSTEM);
				VectorType newVector1(transformSystem.n_row(),0);
				getNewPure(newVector1,pureVectors_.first,ProgramGlobals::SYSTEM,
				           alphaFixed,lrs_.left(),transformSystem,sites.first);
				pureVectors_.first = newVector1;

				const MatrixType& transformEnviron = 
				                        wft_.transform(ProgramGlobals::ENVIRON);
				VectorType newVector2(transformEnviron.n_row(),0);
				getNewPure(newVector2,pureVectors_.second,ProgramGlobals::ENVIRON,
						   betaFixed,lrs_.right(),transformEnviron,sites.second);
				pureVectors_.second = newVector2;
				setFromInfinite(targetVectors_[0],lrs_);
				assert(std::norm(targetVectors_[0])>1e-6);

				systemPrev_.fixed = alphaFixed;
				systemPrev_.permutationInverse = lrs_.left().permutationInverse();
				environPrev_.fixed = betaFixed;
				environPrev_.permutationInverse = lrs_.right().permutationInverse();
			}

			void getFullVector(std::vector<RealType>& v,size_t m,const LeftRightSuperType& lrs) const
			{
				int offset = lrs.super().partition(m);
				int total = lrs.super().partition(m+1) - offset;

				PackIndicesType pack(lrs.left().size());
				v.resize(total);
				assert(PsimagLite::norm(pureVectors_.first)>1e-6);
				assert(PsimagLite::norm(pureVectors_.second)>1e-6);
				for (int i=0;i<total;i++) {
					size_t alpha,beta;
					pack.unpack(alpha,beta,lrs.super().permutation(i+offset));
					v[i] = pureVectors_.first[alpha] * pureVectors_.second[beta];
				}
			}

			void getNewPure(VectorType& newVector,
			                VectorType& oldVector,
			                size_t direction,
			                size_t alphaFixed,
			                const BasisWithOperatorsType& basis,
			                const MatrixType& transform,
			                size_t site)
			{
				if (oldVector.size()==0) setInitialPure(oldVector,site);
				VectorType tmpVector;
				if (transform.n_row()==0) {
					tmpVector = oldVector;
					assert(PsimagLite::norm(tmpVector)>1e-6);
				} else {
					delayedTransform(tmpVector,oldVector,direction,transform,site);
					assert(PsimagLite::norm(tmpVector)>1e-6);
				}
				size_t ns = tmpVector.size();
				size_t nk = model_.hilbertSize(site);
				size_t newSize =  (transform.n_col()==0) ? (ns*ns) : 
							transform.n_col() * nk;
				newVector.resize(newSize);
				for (size_t alpha=0;alpha<newVector.size();alpha++)
					newVector[alpha] = 0;

				for (size_t alpha=0;alpha<ns;alpha++) {
					size_t gamma = (direction==ProgramGlobals::SYSTEM) ?
						basis.permutationInverse(alpha + alphaFixed*ns) :
						basis.permutationInverse(alphaFixed + alpha*nk);
					newVector[gamma] = tmpVector[alpha];
				}
//				for (size_t gamma=0;gamma<newVector.size();gamma++) {
//					newVector[gamma] = 0;
//					for (size_t alpha=0;alpha<ns;alpha++) {
//						size_t gammaPrime = (direction==ProgramGlobals::SYSTEM) ?
//									basis.permutationInverse(alpha + alphaFixed*ns) :
//									basis.permutationInverse(alphaFixed + alpha*nk);

//						if (gamma == gammaPrime)
//							newVector[gamma] += tmpVector[alpha];
//					}
//				}
				std::ostringstream msg2;
				msg2<<"Old size of pure is "<<ns<<" norm="<<PsimagLite::norm(tmpVector);
				progress_.printline(msg2,std::cerr);
				std::ostringstream msg;
				msg<<"New size of pure is "<<newSize<<" norm="<<PsimagLite::norm(newVector);
				progress_.printline(msg,std::cerr);
				assert(PsimagLite::norm(newVector)>1e-6);
			}

			void delayedTransform(VectorType& newVector,
								  VectorType& oldVector,
								  size_t direction,
								  const MatrixType& transform,
								  size_t site)
			{
				assert(oldVector.size()==transform.n_row());

				size_t ne = model_.hilbertSize(site);
				
				const std::vector<size_t>& permutationInverse = (direction==SYSTEM)
				? systemPrev_.permutationInverse : environPrev_.permutationInverse;
				size_t nsPrev = permutationInverse.size()/ne;
				
				newVector.resize(transform.n_col());
				//newVector = oldVector * transform;
				for (size_t gamma=0;gamma<newVector.size();gamma++) {
					newVector[gamma] = 0;
					for (size_t alpha=0;alpha<nsPrev;alpha++) {
						size_t noPermIndex =  (direction==SYSTEM)
								   ? alpha + systemPrev_.fixed*nsPrev
								   : environPrev_.fixed + alpha*ne;
						
						size_t gammaPrime = permutationInverse[noPermIndex];
						
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


				for (size_t gamma=0;gamma<newVector.size();gamma++) {
					newVector[gamma] = 0;
					for (size_t gammaPrime=0;gammaPrime<oldVector.size();gammaPrime++) {
						newVector[gamma] += transform(gammaPrime,gamma) *
									  oldVector[gammaPrime];
					}
				}
			}

			void setInitialPure(VectorType& oldVector,size_t site)
			{
				size_t sitePlusOrMinus = (site==1) ? 0 : site+1;
				size_t alphaFixed = mettsStochastics_.chooseRandomState(sitePlusOrMinus);

				std::ostringstream msg;
				msg<<"New pures for site "<<sitePlusOrMinus<<" is "<<alphaFixed;
				progress_.printline(msg,std::cerr);

				oldVector.resize(model_.hilbertSize(site));
				assert(alphaFixed<oldVector.size());
				for (size_t i=0;i<oldVector.size();i++) {
					oldVector[i] = (i==alphaFixed) ? 1 : 0;
				}
				assert(PsimagLite::norm(oldVector)>1e-6);
			}

			void setFromInfinite(VectorWithOffsetType& phi,const LeftRightSuperType& lrs) const
			{
				phi.populateSectors(lrs.super());
				for (size_t ii=0;ii<phi.sectors();ii++) {
					size_t i0 = phi.sector(ii);
					VectorType v;
					getFullVector(v,i0,lrs);
					RealType tmpNorm = PsimagLite::norm(v);
					if (fabs(tmpNorm-1.0)<1e-6) {
						size_t j = lrs.super().qn(lrs.super().partition(i0));
						std::vector<size_t> qns = BasisType::decodeQuantumNumber(j);
						std::cerr<<"setFromInfinite: qns= ";
						for (size_t k=0;k<qns.size();k++) std::cerr<<qns[k]<<" ";
						std::cerr<<"\n";
					}
					phi.setDataInSector(v,i0);
				}
				phi.collapseSectors();
				assert(std::norm(phi)>1e-6);
			}

			// in situ computation:
			void cocoon(size_t direction,const PairType& sites,bool corner)
			{
				std::cerr<<"-------------&*&*&* In-situ measurements start\n";
				//test(psi_,psi_,direction,"<PSI|A|PSI>",sites);
				
				for (size_t j=0;j<targetVectors_.size();j++) {
					std::string s = "<P"+ttos(j)+"|A|P"+ttos(j)+">";
					size_t site = test(targetVectors_[j],targetVectors_[j],direction,s,sites,corner);
					if (stage_==COLLAPSE && j==0) sitesCollapsed_.push_back(site);
				}
				std::cerr<<"-------------&*&*&* In-situ measurements end\n";
			}

			void checkOrder(size_t i) const
			{
				if (i==0) return;
				for (size_t j=0;j<i;j++) {
					if (stage_ == DISABLED) {
						std::string s ="TST:: Seeing tst site "+ttos(mettsStruct_.sites[i]);
						s =s + " before having seen";
						s = s + " site "+ttos(j);
						s = s +". Please order your tst sites in order of appearance.\n";
						throw std::runtime_error(s);
					}
				}
			}

			bool allStages(size_t x) const
			{
				if (stage_!=x) return false;
				return true;
			}

			bool noStageIs(size_t x) const
			{
				//for (size_t i=0;i<stage_.size();i++)
				if (stage_==x) return false;
				return true;
			}

			std::string getStage() const
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
			                     size_t systemOrEnviron)
			{
				const VectorWithOffsetType& phi = targetVectors_[startEnd.first];
				std::ostringstream msg;
				msg<<" vector number "<<startEnd.first<<" has norm ";
				msg<<std::norm(phi);
				progress_.printline(msg,std::cout);
				if (std::norm(phi)<1e-6) setFromInfinite(targetVectors_[startEnd.first],lrs_);
				calcTimeVectorsTest(startEnd,Eg,phi,systemOrEnviron);
			}

			void calcTimeVectorsTest(const PairType& startEnd,
			                            RealType Eg,
			                            const VectorWithOffsetType& phi,
			                            size_t systemOrEnviron)
			{
				std::vector<MatrixType> V(phi.sectors());
				std::vector<MatrixType> T(phi.sectors());

				std::vector<size_t> steps(phi.sectors());

				triDiag(phi,T,V,steps);

				std::vector<std::vector<RealType> > eigs(phi.sectors());
						
				for (size_t ii=0;ii<phi.sectors();ii++) 
					PsimagLite::diag(T[ii],eigs[ii],'V');
				
				calcTargetVectors(startEnd,phi,T,V,Eg,eigs,steps,systemOrEnviron);
			}

			void calcTargetVectors(const PairType& startEnd,
			                       const VectorWithOffsetType& phi,
			                       const std::vector<MatrixType>& T,
			                       const std::vector<MatrixType>& V,
			                       RealType Eg,
			                       const std::vector<VectorType>& eigs,
					       std::vector<size_t> steps,
			                       size_t systemOrEnviron)
			{
				for (size_t i=startEnd.first+1;i<startEnd.second;i++) {
					VectorWithOffsetType v;
					// Only time differences here (i.e. betas_[i] not betas_[i]+currentBeta_)
					calcTargetVector(v,phi,T,V,Eg,eigs,i,steps);
					RealType x = 1.0/std::norm(v);
					targetVectors_[i]= x* v;
				}
			}

			void calcTargetVector(VectorWithOffsetType& v,
			                      const VectorWithOffsetType& phi,
			                      const std::vector<MatrixType>& T,
			                      const std::vector<MatrixType>& V,
			                      RealType Eg,
			                      const std::vector<VectorType>& eigs,
	    	                      size_t timeIndex,
			                      std::vector<size_t> steps)
			{
				v = phi;
				for (size_t ii=0;ii<phi.sectors();ii++) {
					size_t i0 = phi.sector(ii);
					VectorType r;
					calcTargetVector(r,phi,T[ii],V[ii],Eg,eigs[ii],timeIndex,steps[ii],i0);
					v.setDataInSector(r,i0);
				}
			}

			void calcTargetVector(VectorType& r,
      		                      const VectorWithOffsetType& phi,
			                      const MatrixType& T,
			                      const MatrixType& V,
			                      RealType Eg,
			                      const VectorType& eigs,
			                      size_t timeIndex,
			                      size_t steps,
			                      size_t i0)
			{
				size_t n2 = steps;
				size_t n = V.n_row();
				assert(T.n_col()==T.n_row());
				assert(V.n_col()==T.n_col());

				RealType rone = 1.0;
				RealType rzero = 0.0;
				
				VectorType tmp(n2);
				r.resize(n2);
				calcR(r,T,V,phi,Eg,eigs,timeIndex,steps,i0);
				psimag::BLAS::GEMV('N', n2, n2, rone, &(T(0,0)), n2, &(r[0]), 1, rzero, &(tmp[0]), 1 );
				r.resize(n);
				psimag::BLAS::GEMV('N', n,  n2, rone, &(V(0,0)), n, &(tmp[0]),1, rzero, &(r[0]),   1 );
			}

			void calcR(VectorType& r,
				   const MatrixType& T,
			           const MatrixType& V,
				   const VectorWithOffsetType& phi,
				   RealType Eg,
				   const VectorType& eigs,
				   size_t timeIndex,
			           size_t n2,
			           size_t i0)
			{
				for (size_t k=0;k<n2;k++) {
					RealType sum = 0.0;
					for (size_t kprime=0;kprime<n2;kprime++) {
						RealType tmpV = calcVTimesPhi(kprime,V,phi,i0);
						sum += std::conj(T(kprime,k))*tmpV;
					}
					RealType tmp = (eigs[k])*betas_[timeIndex];
					r[k] = sum * exp(-tmp);
				}
			}

			RealType calcVTimesPhi(size_t kprime,
			                      const MatrixType& V,
			                      const VectorWithOffsetType& phi,
			                      size_t i0)
			{
				RealType ret = 0;
				size_t total = phi.effectiveSize(i0);
				
				for (size_t j=0;j<total;j++)
					ret += std::conj(V(j,kprime))*phi.fastAccess(i0,j);
				return ret;
			}

			void triDiag(const VectorWithOffsetType& phi,
			             std::vector<MatrixType>& T,
	 		             std::vector<MatrixType>& V,
			             std::vector<size_t>& steps)
			{
				for (size_t ii=0;ii<phi.sectors();ii++) {
					size_t i = phi.sector(ii);
					steps[ii] = triDiag(phi,T[ii],V[ii],i);
				}
			}

			size_t triDiag(const VectorWithOffsetType& phi,
			               MatrixType& T,
			               MatrixType& V,
			               size_t i0)
			{
				size_t p = lrs_.super().findPartitionNumber(phi.offset(i0));
				typename ModelType::ModelHelperType modelHelper(p,lrs_);
				 		//,useReflection_);
				typename LanczosSolverType::LanczosMatrixType lanczosHelper(&model_,&modelHelper);

				ParametersForSolverType params;
				params.steps = model_.params().lanczosSteps;
				params.tolerance = model_.params().lanczosEps;
				params.stepsForEnergyConvergence =ProgramGlobals::MaxLanczosSteps;

				LanczosSolverType lanczosSolver(lanczosHelper,params,&V);

				TridiagonalMatrixType ab;
				size_t total = phi.effectiveSize(i0);
				TargetVectorType phi2(total);
				phi.extract(phi2,i0);
				RealType x = PsimagLite::norm(phi2);
				assert(x>1e-6);
				std::cerr<<"norm of phi2="<<x<<"\n";
				lanczosSolver.decomposition(phi2,ab);
				lanczosSolver.buildDenseMatrix(T,ab);
				return lanczosSolver.steps();
			}

			RealType setOneInterval(const RealType& factor,
			                        const PairType& startEnd,
			                        const RealType& tau)
			{
				RealType sum = 0;
				for (size_t i=startEnd.first;i<startEnd.second;i++) {
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

			void findElectronsOfOneSite(std::vector<size_t>& electrons,size_t site) const
			{
				std::vector<size_t> block(1,site);
				typename ModelType::HilbertBasisType basis;
				std::vector<size_t> quantumNumbs;
				model_.setNaturalBasis(basis,quantumNumbs,block);
				model_.findElectrons(electrons,basis,site);
			}

			size_t test(const VectorWithOffsetType& src1,
			          const VectorWithOffsetType& src2,
			          size_t systemOrEnviron,
			          const std::string& label,
					  const PairType& sites,
					  bool corner) const
			{
				VectorWithOffsetType dest;
				OperatorType A = getObservableToTest(model_.params().model);
//				typename ModelType::HilbertBasisType basis;
//				std::vector<size_t> quantumNumbs;
//				assert(sites.first==sites.second);
				std::vector<size_t> electrons;
				findElectronsOfOneSite(electrons,sites.first);
				FermionSign fs(lrs_.left(),electrons);
				applyOpLocal_(dest,src1,A,fs,systemOrEnviron,corner);

				RealType sum = 0;
				for (size_t ii=0;ii<dest.sectors();ii++) {
					size_t i = dest.sector(ii);
					size_t offset1 = dest.offset(i);
					for (size_t jj=0;jj<src2.sectors();jj++) {
						size_t j = src2.sector(jj);
						size_t offset2 = src2.offset(j);
						if (i!=j) continue; //throw std::runtime_error("Not same sector\n");
						for (size_t k=0;k<dest.effectiveSize(i);k++) 
							sum+= dest[k+offset1] * std::conj(src2[k+offset2]);
					}
				}
				RealType nor = std::norm(src1);
				size_t site = sites.first;
				if (corner) {
					if (site==1) {
						site=0;
					} else {
						site=model_.geometry().numberOfSites()-1;
					}
				}

				std::cerr<<site<<" "<<sum<<" "<<" "<<currentBeta_;
				std::cerr<<" "<<label<<" "<<nor<<" "<<std::norm(src2);
				std::cerr<<" "<<std::norm(dest)<<"    "<<sum/(nor*nor)<<"\n";

				return site;
			}

			OperatorType getObservableToTest(const std::string& modelName) const
			{
				OperatorType A;
				size_t site = 0; // sites.first; <-- site-dependent Hilbert space not supported by METTS

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
				std::string s(__FILE__);
				s += " " + ttos(__LINE__) + "\n";
				s += "Model " + modelName + " not supported by MettsTargetting\n";
				throw std::runtime_error(s.c_str());
			}

			void printEnergies() const
			{
				for (size_t i=0;i<targetVectors_.size();i++)
					printEnergies(targetVectors_[i],i);
			}

			void printEnergies(const VectorWithOffsetType& phi,size_t whatTarget) const
			{
				for (size_t ii=0;ii<phi.sectors();ii++) {
					size_t i = phi.sector(ii);
					printEnergies(phi,whatTarget,i);
				}
			}

			void printEnergies(const VectorWithOffsetType& phi,size_t whatTarget, size_t i0) const
			{
				size_t p = lrs_.super().findPartitionNumber(phi.offset(i0));
				typename ModelType::ModelHelperType modelHelper(p,lrs_);
						//,useReflection_);
				typename LanczosSolverType::LanczosMatrixType lanczosHelper(&model_,&modelHelper);


				size_t total = phi.effectiveSize(i0);
				TargetVectorType phi2(total);
				phi.extract(phi2,i0);
				TargetVectorType x(total);
				lanczosHelper.matrixVectorProduct(x,phi2);
				std::ostringstream msg;
				msg<<"Hamiltonian average at beta="<<currentBeta_<<" for target="<<whatTarget;
				msg<<" sector="<<i0<<" <phi(t)|H|phi(t)>="<<(phi2*x)<<" <phi(t)|phi(t)>="<<(phi2*phi2);
				progress_.printline(msg,std::cout);
			}

			bool isAtBorder(size_t direction,size_t site) const
			{
				if (direction==EXPAND_SYSTEM && site+2==model_.geometry().numberOfSites())
					return true;
				if (direction!=EXPAND_SYSTEM && site==1)
					return true;
				return false;
			}

			bool allSitesCollapsed() const
			{
				size_t n = model_.geometry().numberOfSites();
				for (size_t i=0;i<n;i++) {
					bool seen = (std::find(sitesCollapsed_.begin(),sitesCollapsed_.end(),i) != sitesCollapsed_.end());
					if (!seen) return false;
				}
				return true;
			}

			size_t stage_;
			const LeftRightSuperType& lrs_;
			const ModelType& model_;
			const TargettingParamsType& mettsStruct_;
			const WaveFunctionTransfType& wft_;
			const size_t& quantumSector_;
			PsimagLite::ProgressIndicator progress_;
			RealType currentBeta_;
			std::vector<RealType> betas_,weight_;
			std::vector<VectorWithOffsetType> targetVectors_;
			RealType gsWeight_;
			ApplyOperatorType applyOpLocal_;
			MettsStochasticsType mettsStochastics_;
			MettsCollapseType mettsCollapse_;
			size_t timesWithoutAdvancement_;
			size_t prevDirection_;
			MettsPrev systemPrev_;
			MettsPrev environPrev_;
			std::pair<VectorType,VectorType> pureVectors_;
			std::vector<size_t> sitesCollapsed_;
	};     //class MettsTargetting

	template<
		template<typename,typename,typename> class LanczosSolverTemplate,
			template<typename,typename> class InternalProductTemplate,
 		template<typename,typename> class WaveFunctionTransfTemplate,
			typename ModelType_,
 		typename ConcurrencyType_,
			typename IoType_,
   			template<typename> class VectorWithOffsetTemplate>
	std::ostream& operator<<(std::ostream& os,
			const MettsTargetting<LanczosSolverTemplate,
			InternalProductTemplate,
			WaveFunctionTransfTemplate,ModelType_,ConcurrencyType_,IoType_,
			VectorWithOffsetTemplate>& tst)
	{
		tst.print(os);
		return os;
	}
} // namespace Dmrg

#endif //METTS_TARGETTING_H
