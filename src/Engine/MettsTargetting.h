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
			typedef std::vector<RealType> VectorType;
			typedef PsimagLite::ParametersForSolver<RealType> ParametersForSolverType;
			typedef LanczosSolverTemplate<ParametersForSolverType,InternalProductType,VectorType>
			    LanczosSolverType;
			typedef typename LanczosSolverType::TridiagonalMatrixType TridiagonalMatrixType;
			typedef typename BasisWithOperatorsType::OperatorType OperatorType;
			typedef typename BasisWithOperatorsType::BasisType BasisType;
			typedef MettsParams<ModelType> TargettingParamsType;
			typedef typename BasisType::BlockType BlockType;
			typedef VectorWithOffsetTemplate<RealType> VectorWithOffsetType;
			typedef WaveFunctionTransfTemplate<LeftRightSuperType,VectorWithOffsetType> WaveFunctionTransfType;
			typedef VectorType TargetVectorType;
			typedef PsimagLite::Matrix<RealType> MatrixType;
			typedef BlockMatrix<RealType,MatrixType> BlockMatrixType;
			typedef ApplyOperatorLocal<LeftRightSuperType,VectorWithOffsetType,TargetVectorType> ApplyOperatorType;
			typedef MettsSerializer<RealType,VectorWithOffsetType> MettsSerializerType;
			typedef MettsStochastics<ModelType> MettsStochasticsType;
			typedef typename MettsStochasticsType::PairType PairType;
			typedef MettsCollapse<VectorWithOffsetType,MettsStochasticsType> MettsCollapseType;
			typedef typename MettsCollapseType::PackIndicesType PackIndicesType;
			
			enum {DISABLED,WFT_NOADVANCE,WFT_ADVANCE};
			enum {EXPAND_ENVIRON=WaveFunctionTransfType::EXPAND_ENVIRON,
			EXPAND_SYSTEM=WaveFunctionTransfType::EXPAND_SYSTEM,
			INFINITE=WaveFunctionTransfType::INFINITE};
			const static size_t SYSTEM = ProgramGlobals::SYSTEM;

			// static size_t const PRODUCT = TargettingParamsType::PRODUCT;
			// static size_t const SUM = TargettingParamsType::SUM;
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
			  mettsCollapse_(mettsStochastics_,lrs_),
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
				return (allStages(DISABLED)) ? 1 : targetVectors_.size();
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

				// Advance or wft  collapsed vector
				evolve(n1,n1,n1-1,Eg,direction,sites,loopNumber);

				// compute imag. time evolution:
				calcTimeVectors(PairType(0,n1),Eg,direction);

				// collapse
				bool hasCollapsed = mettsCollapse_(targetVectors_[n1],
				                       targetVectors_[0],sites.first,direction);
				
				if (hasCollapsed) {
					std::string s = "  COLLAPSEHERE  ";
					test(targetVectors_[n1],targetVectors_[n1],direction,s,sites);
				}

				// in-situ measurement
				cocoon(direction,sites); 
			}

			void load(const std::string& f)
			{
				throw std::runtime_error("Metts: load() unimplemented\n");
// 				for (size_t i=0;i<stage_.size();i++) stage_[i] = WFT_NOADVANCE;
// 
// 				typename IoType::In io(f);
// 
// 				TimeSerializerType ts(io,IoType::In::LAST_INSTANCE);
// 				for (size_t i=0;i<targetVectors_.size();i++) targetVectors_[i] = ts.vector(i);
// 				currentBeta_ = ts.time();
// 
// 				psi_.load(io,"PSI");
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
// 				std::ostringstream msg;
// 				msg<<"WARNING: initial guess: Needs work";
// 				progress_.printline(msg,std::cout);
//  				wft_.setInitialVector(v,psi_,lrs_);
				
// 				bool b = allStages(WFT_ADVANCE) || allStages(WFT_NOADVANCE);
// 				if (!b) return;
// 				std::vector<VectorWithOffsetType> vv(targetVectors_.size());
// 				for (size_t i=0;i<targetVectors_.size();i++) {
// 					wft_.setInitialVector(vv[i],
// 						targetVectors_[i],lrs_);
// 					if (norm(vv[i])<1e-6) continue;
// 					VectorWithOffsetType w= weight_[i]*vv[i];
// 					v += w;
// 				}
			}

			template<typename IoOutputType>
			void save(const std::vector<size_t>& block,IoOutputType& io) const
			{
				std::ostringstream msg;
				msg<<"WARNING: save(...) unimplemented";
// 				msg<<"Saving state...";
 				progress_.printline(msg,std::cout);
// 
// 				TimeSerializerType ts(currentBeta_,block[0],targetVectors_);
// 				ts.save(io);
// 				psi_.save(io,"PSI");
			}

			RealType time() const { return 0; }

			void updateOnSiteForTimeDep(BasisWithOperatorsType& basisWithOps) const
			{
				// nothing to do here
			}

//			void truncate(BasisWithOperatorsType& pS,BasisWithOperatorsType& pE)
//			{
//				const MatrixType& transformSystem =  wft_.transform(ProgramGlobals::SYSTEM);
//				VectorType newVector1 = pureVectors_.first;
//				pureVectors_.first =  newVector1 * transformSystem ;

//				const MatrixType& transformEnviron =
//						wft_.transform(ProgramGlobals::ENVIRON);
//				VectorType newVector2 = pureVectors_.second;
//				pureVectors_.second = newVector2 * transformEnviron;

//				BasisType super("super");
//				super.setToProduct(pS,pE);
//				LeftRightSuperType lrs(pS,pE,super);
//				setFromInfinite(targetVectors_[0],lrs);

//				assert(std::norm(targetVectors_[0])>1e-6);

//				std::ostringstream msg;
//				msg<<"Truncating, targetVectors_[0].size="<<targetVectors_[0].size();
//				progress_.printline(msg,std::cerr);
//			}

		private:

			void evolve(size_t index,
			            size_t start,
			            size_t indexAdvance,
			            RealType Eg,
			            size_t direction,
			            std::pair<size_t,size_t> sites,
			            size_t loopNumber)
			{
				if (index==0 && start==0) advanceCounterAndComputeStage();

				std::ostringstream msg;
				msg<<"Evolving, stage="<<getStage()<<" loopNumber="<<loopNumber;
				msg<<" Eg="<<Eg;
				progress_.printline(msg,std::cout);
				assert(sites.first==sites.second);
				size_t nk = model_.hilbertSize(sites.first);
				advanceOrWft(index,indexAdvance,direction,nk,sites);
			}

			void advanceCounterAndComputeStage()
			{
				stage_=WFT_NOADVANCE;

				if (timesWithoutAdvancement_ >= mettsStruct_.advanceEach) {
					stage_ = WFT_ADVANCE;
					currentBeta_ += mettsStruct_.tau;
					timesWithoutAdvancement_=0;
				} else {
					if (stage_==WFT_NOADVANCE) timesWithoutAdvancement_++;
				}

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
// 				size_t indexAdvance = betas_.size()-1; // FIXME 
// 				size_t indexNoAdvance = 0;

				if (stage_== WFT_NOADVANCE || stage_== WFT_ADVANCE) {
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

					examineTarget(targetVectors_[advance],advance);
					//populateCorrectSector(phiNew,lrs_);
					phiNew.populateSectors(lrs_.super());
					// OK, now that we got the partition number right, let's wft:
					wft_.setInitialVector(phiNew,targetVectors_[advance],lrs_,nk);
					phiNew.collapseSectors();
					assert(std::norm(phiNew)>1e-6);
					targetVectors_[index] = phiNew;
					//examineTarget(targetVectors_[index],index);
				} else {
					assert(false);
				}
			}

			void updateStochastics(const PairType& sites)
			{
				// only way of getting the quantum number where the
				//  the pure resides
				size_t m = getPartition();
				size_t qn = lrs_.super().qn(lrs_.super().partition(m));
				
				mettsStochastics_.update(qn,sites);
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
				msg<<"New pures for site"<<sites;
				progress_.printline(msg,std::cerr);

				const MatrixType& transformSystem =  wft_.transform(ProgramGlobals::SYSTEM);
				VectorType newVector1(transformSystem.n_row(),0);
				getNewPure(newVector1,pureVectors_.first,ProgramGlobals::SYSTEM,
				           alphaFixed,lrs_.left(),transformSystem,sites.first);
// 				systemPrev_.ns = pureVectors_.first.size();
				pureVectors_.first = newVector1;

				const MatrixType& transformEnviron = 
				                        wft_.transform(ProgramGlobals::ENVIRON);
				VectorType newVector2(transformEnviron.n_row(),0);
				getNewPure(newVector2,pureVectors_.second,ProgramGlobals::ENVIRON,
						   betaFixed,lrs_.right(),transformEnviron,sites.second);
// 				environPrev_.ns = pureVectors_.second.size();
				pureVectors_.second = newVector2;
				setFromInfinite(targetVectors_[0],lrs_);
				//examineTarget(targetVectors_[0],0);
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
				} else {
					delayedTransform(tmpVector,oldVector,direction,transform,site);
				}
				size_t ns = tmpVector.size();
				size_t nk = model_.hilbertSize(site);
				size_t newSize =  (transform.n_col()==0) ? (ns*ns) : 
							transform.n_col() * nk;
				newVector.resize(newSize);
				//for (size_t alpha=0;alpha<newVector.size();alpha++) newVector[alpha] = 0;

//				for (size_t alpha=0;alpha<ns;alpha++) {
//					size_t gamma = (direction==ProgramGlobals::SYSTEM) ?
//					    basis.permutationInverse(alpha + alphaFixed*ns) :
//					    basis.permutationInverse(alphaFixed + alpha*nk);
//					newVector[gamma] = tmpVector[alpha];
//				}
				for (size_t gamma=0;gamma<newVector.size();gamma++) {
					newVector[gamma] = 0;
					for (size_t alpha=0;alpha<ns;alpha++) {
						size_t gammaPrime = (direction==ProgramGlobals::SYSTEM) ?
									basis.permutationInverse(alpha + alphaFixed*ns) :
									basis.permutationInverse(alphaFixed + alpha*nk);

						if (gamma == gammaPrime)
							newVector[gamma] += tmpVector[alpha];
					}
				}
				std::ostringstream msg;
				msg<<"New size of pure is "<<newSize<<" norm="<<PsimagLite::norm(newVector);
				progress_.printline(msg,std::cerr);
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

			void setInitialPure(VectorType& oldVector,size_t site)
			{
				size_t sitePlusOrMinus = (site==1) ? 0 : site+1;
				size_t alphaFixed = mettsStochastics_.chooseRandomState(sitePlusOrMinus);
				oldVector.resize(model_.hilbertSize(site));
				for (size_t i=0;i<oldVector.size();i++) {
					oldVector[i] = (i==alphaFixed) ? 1 : 0;
				}
			}

//			void populateCorrectSector(VectorWithOffsetType& phi,const LeftRightSuperType& lrs) const
//			{
//				size_t total = lrs.super().partition()-1;
//				size_t m = getPartition();
//				std::vector<VectorType> vv;
//				for (size_t i=0;i<total;i++) {
//					size_t bs = lrs.super().partition(i+1)-lrs.super().partition(i);
//					if (i!=m) bs=0;
//					VectorType vone(bs);
//					vv.push_back(vone);
//				}
				
//				phi.set(vv,lrs.super());
//				assert(phi.sectors()==1);
//			}

			void setFromInfinite(VectorWithOffsetType& phi,const LeftRightSuperType& lrs) const
			{
				//populateCorrectSector(phi,lrs);
				phi.populateSectors(lrs.super());
				std::cerr<<"NUUUUUUMBBBBBBBERRRRRRRRR OF SEEEEECTTTTTTTTORRRRRRRRS="<<phi.sectors()<<"\n";
				for (size_t ii=0;ii<phi.sectors();ii++) {
					size_t i0 = phi.sector(ii);
					VectorType v;
					getFullVector(v,i0,lrs);
					phi.setDataInSector(v,i0);
				}
				phi.collapseSectors();
				assert(std::norm(phi)>1e-6);
			}

			// in situ computation:
			void cocoon(size_t direction,const PairType& sites) const
			{
				std::cerr<<"-------------&*&*&* In-situ measurements start\n";
				//test(psi_,psi_,direction,"<PSI|A|PSI>",sites);
				
				for (size_t j=0;j<targetVectors_.size();j++) {
					std::string s = "<P"+ttos(j)+"|A|P"+ttos(j)+">";
					test(targetVectors_[j],targetVectors_[j],direction,s,sites);
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
// 					case OPERATOR:
// 						return "Applying operator for the first time";
// 						break; 
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

				std::vector<MatrixType> V(phi.sectors());
				std::vector<MatrixType> T(phi.sectors());

				std::vector<size_t> steps(phi.sectors());

				triDiag(phi,T,V,steps);

				std::vector<std::vector<RealType> > eigs(phi.sectors());
						
				for (size_t ii=0;ii<phi.sectors();ii++) 
					PsimagLite::diag(T[ii],eigs[ii],'V');
				
				calcTargetVectors(startEnd,T,V,Eg,eigs,steps,systemOrEnviron);
			}

			void calcTargetVectors(const PairType& startEnd,
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
					calcTargetVector(v,targetVectors_[startEnd.first],T,V,Eg,eigs,i,steps);
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
					RealType tmp = (eigs[k]-Eg)*betas_[timeIndex];
					r[k] = sum * exp(tmp);
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
			
				RealType eps= 0.01*ProgramGlobals::LanczosTolerance;
				size_t iter= ProgramGlobals::LanczosSteps;

				ParametersForSolverType params;
				params.steps = iter;
				params.tolerance = eps;
				params.stepsForEnergyConvergence =ProgramGlobals::MaxLanczosSteps;

				LanczosSolverType lanczosSolver(lanczosHelper,params,&V);

				TridiagonalMatrixType ab;
				size_t total = phi.effectiveSize(i0);
				TargetVectorType phi2(total);
				phi.extract(phi2,i0);
				/* std::ostringstream msg;
				msg<<"Calling tridiagonalDecomposition...\n";
				progress_.printline(msg,std::cerr);*/
				RealType x = PsimagLite::norm(phi2);
				assert(x>1e-6);
				std::cerr<<"norm of phi2="<<x<<"\n";
				lanczosSolver.decomposition(phi2,ab);
				lanczosSolver.buildDenseMatrix(T,ab);
				//check1(V,phi2);
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

			//! This check is invalid if there are more than one sector
			void check1(const MatrixType& V,const TargetVectorType& phi2)
			{
				assert(V.n_col()<=V.n_row());
				TargetVectorType r(V.n_col());
				for (size_t k=0;k<V.n_col();k++) {
					r[k] = 0.0;
					for (size_t j=0;j<V.n_row();j++) 
						r[k] += conj(V(j,k))*phi2[j];
					// is r(k) == \delta(k,0)
					if (k==0 && std::norm(r[k]-1.0)>1e-5) 
						std::cerr<<"WARNING: r[0]="<<r[0]<<" != 1\n";
					if (k>0 && std::norm(r[k])>1e-5) 
						std::cerr<<"WARNING: r["<<k<<"]="<<r[k]<<" !=0\n";
				}
			}

// 			void guessPhiSectors(VectorWithOffsetType& phi,size_t i,size_t systemOrEnviron)
// 			{
// 				FermionSign fs(lrs_.left(),mettsStruct_.electrons);
// 				if (allStages(WFT_NOADVANCE)) {
// 					VectorWithOffsetType tmpVector = psi_;
// 					for (size_t j=0;j<mettsStruct_.aOperators.size();j++) {
// 						applyOpLocal_(phi,tmpVector,mettsStruct_.aOperators[j],fs,
// 							systemOrEnviron);
// 						tmpVector = phi;
// 					}
// 					return;
// 				}
// 				applyOpLocal_(phi,psi_,mettsStruct_.aOperators[i],fs,
// 								systemOrEnviron);
// 			}

			void examineTarget(const VectorWithOffsetType& phi,size_t index) const
			{
				std::cerr<<"HEEEEEEEREEEEEEEE size="<<phi.size()<<" index="<<index<<" sectors="<<phi.sectors()<<"\n";
				for (size_t ii=0;ii<phi.sectors();ii++) {
					size_t i = phi.sector(ii);
					VectorType phi2;
					phi.extract(phi2,i);
					std::cerr<<"NOOOOOOORMMMMMMMM of "<<i<<" is "<<PsimagLite::norm(phi2)<<" index="<<index;
					size_t qn = lrs_.super().qn(phi.offset(i));
					std::cerr<<" sectorSize="<<phi2.size()<<" QN="<<qn<<"\n";
					PsimagLite::IoSimple::Out io("test1.txt",0);
					io.printVector(phi2,"phi2");
				}
			}

			void zeroOutVectors()
			{
				for (size_t i=0;i<targetVectors_.size();i++) 
					targetVectors_[i].resize(lrs_.super().size());
			}

			void findElectronsOfOneSite(std::vector<size_t>& electrons,size_t site) const
			{
				std::vector<size_t> block(1,site);
				typename ModelType::HilbertBasisType basis;
				std::vector<size_t> quantumNumbs;
				model_.setNaturalBasis(basis,quantumNumbs,block);
				model_.findElectrons(electrons,basis,site);
			}

			void test(const VectorWithOffsetType& src1,
			          const VectorWithOffsetType& src2,
			          size_t systemOrEnviron,
			          const std::string& label,
			          const PairType& sites) const
			{
				VectorWithOffsetType dest;
				OperatorType A;
				size_t site = 0; // sites.first; <-- site-dependent Hilbert space not supported by METTS
				PsimagLite::CrsMatrix<RealType> tmpC(model_.naturalOperator("nup",site,0));
				A.data = tmpC;
//				PsimagLite::CrsMatrix<RealType> tmpCt;
//				transposeConjugate(tmpCt,tmpC);
//				multiply(A.data,tmpCt,tmpC);
				A.fermionSign = 1;
				//A.data = tmpC;
				typename ModelType::HilbertBasisType basis;
				std::vector<size_t> quantumNumbs;
//				assert(sites.first==sites.second);
				std::vector<size_t> electrons;
				findElectronsOfOneSite(electrons,sites.first);
				FermionSign fs(lrs_.left(),electrons);
				applyOpLocal_(dest,src1,A,fs,systemOrEnviron);

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
				std::cerr<<sites.first<<" "<<sum<<" "<<" "<<currentBeta_;
				std::cerr<<" "<<label<<" "<<nor<<" "<<std::norm(src2);
				std::cerr<<" "<<std::norm(dest)<<"    "<<sum/(nor*nor)<<"\n";
			}

			size_t stage_;
			//VectorWithOffsetType psi_;
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
			//typename IoType::Out io_;
			ApplyOperatorType applyOpLocal_;
			MettsStochasticsType mettsStochastics_;
			MettsCollapseType mettsCollapse_;
			size_t timesWithoutAdvancement_;
			size_t prevDirection_;
			MettsPrev systemPrev_;
			MettsPrev environPrev_;
			std::pair<VectorType,VectorType> pureVectors_;
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
