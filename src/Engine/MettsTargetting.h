
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
#include "PackIndices.h"
#include "MettsStochastics.h"

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
			typedef PsimagLite::PackIndices PackIndicesType;

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
			typedef LanczosSolverTemplate<RealType,InternalProductType,VectorType>
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
			
			enum {DISABLED,WFT_NOADVANCE,WFT_ADVANCE};
			enum {EXPAND_ENVIRON=WaveFunctionTransfType::EXPAND_ENVIRON,
			EXPAND_SYSTEM=WaveFunctionTransfType::EXPAND_SYSTEM,
			INFINITE=WaveFunctionTransfType::INFINITE};

			static size_t const PRODUCT = TargettingParamsType::PRODUCT;
			static size_t const SUM = TargettingParamsType::SUM;
			static const size_t parallelRank_ = 0; // Metts needs to support concurrency FIXME

			MettsTargetting(const LeftRightSuperType& lrs,
	 		                const ModelType& model,
			                const TargettingParamsType& mettsStruct,
			                const WaveFunctionTransfType& wft)
			: stage_(DISABLED),
			  lrs_(lrs),
			  model_(model),
			  mettsStruct_(mettsStruct),
			  wft_(wft),
			  progress_("MettsTargetting",0),
			  currentBeta_(0),
			  applyOpLocal_(lrs),
			  hilbertSizePerSite_(model_.hilbertSize()),
			  mettsStochastics_(model),
			  timesWithoutAdvancement_(0)
			{
				if (!wft.isEnabled()) throw std::runtime_error(" MettsTargetting "
							"needs an enabled wft\n");
				
				RealType tau =mettsStruct_.tau/(mettsStruct_.timeSteps-1);
				size_t n1 = size_t(mettsStruct_.timeSteps/2);
				if (mettsStruct_.timeSteps & 1) n1++;
				size_t n = mettsStruct_.timeSteps + n1;
				
				betas_.resize(n);
				weight_.resize(n);
				targetVectors_.resize(n);
				
				gsWeight_= 0;
				RealType factor = (1.0 - gsWeight_)/(n+4);
				RealType sum = setOneInterval(factor,PairType(0,n1),tau*0.5);
				sum += setOneInterval(factor,PairType(n1,n),tau);
				
				
				sum += gsWeight_;
				//for (size_t i=0;i<weight_.size();i++) sum += weight_[i];
				if (fabs(sum-1.0)>1e-5)
					throw std::runtime_error("Weights don't amount to one\n");
				//printHeader();
			}

			RealType weight(size_t i) const
			{
				if (allStages(DISABLED)) 
					return 0.5;
				return weight_[i];
			}

			RealType gsWeight() const
			{
				return (allStages(DISABLED)) ?  0.5 : gsWeight_;
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
				psi_.set(v,someBasis);
			}

			const RealType& operator[](size_t i) const { return psi_[i]; }
			
			RealType& operator[](size_t i) { return psi_[i]; }

			const VectorWithOffsetType& gs() const { return psi_; }

			bool includeGroundStage() const {return true; }

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
				if (block1.size()!=1) throw std::runtime_error(
					"MettsTargetting::evolve(...): blocks of size != 1 are unsupported (sorry)\n");
				
				PairType sites(block1[0],block2[0]);
				if (direction==INFINITE) {
					getNewPures(sites);
					return;
				}
				size_t n1 = size_t(mettsStruct_.timeSteps/2);
				if (mettsStruct_.timeSteps & 1) n1++;
				size_t n = mettsStruct_.timeSteps + n1;

				// Advance or wft each target vector for beta/2
				for (size_t i=0;i<n1;i++) {
					evolve(i,0,Eg,direction,sites,loopNumber);
				}
				
				// Advance or wft each target vector for beta
				for (size_t i=n1;i<n;i++) {
					evolve(i,n1,Eg,direction,sites,loopNumber);
				}
				
				calcTimeVectors(PairType(0,n1),Eg,direction);
				calcTimeVectors(PairType(n1,n),Eg,direction);
				
				cocoon(direction,sites); // in-situ
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

			void initialGuess(VectorWithOffsetType& v) const
			{
				std::ostringstream msg;
				msg<<"WARNING: initial guess: Needs work";
				progress_.printline(msg,std::cout);
 				wft_.setInitialVector(v,psi_,lrs_);
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

		private:

			void evolve(size_t index,
						size_t start,
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
				advanceOrWft(index,start,direction);
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
			                  size_t start,
			                  size_t systemOrEnviron)
			{
				if (targetVectors_[index].size()==0) return;

// 				size_t indexAdvance = betas_.size()-1; // FIXME 
// 				size_t indexNoAdvance = 0;
				if (stage_== WFT_ADVANCE) 
					throw std::runtime_error("Advance unimplemented for Metts\n");
				if (stage_== WFT_NOADVANCE || stage_== WFT_ADVANCE) {
// 					size_t advance = indexNoAdvance;
// 					if (stage_[i] == WFT_ADVANCE) advance = indexAdvance;
					std::ostringstream msg;
					msg<<"I'm calling the WFT now";
					progress_.printline(msg,std::cout);

					VectorWithOffsetType phiNew = psi_; // same sectors as g.s.

					// OK, now that we got the partition number right, let's wft:
					wft_.setInitialVector(phiNew,targetVectors_[index],lrs_);
					phiNew.collapseSectors(); 
					targetVectors_[index] = phiNew;
				} else {
					throw std::runtime_error("It's 5 am, do you know what line "
					" your code is exec-ing?\n");
				}
			}

			void getNewPures(const PairType& sites)
			{
				size_t m = psi_.sector(0);
				
				// N.B.: it's really BEFORE_TRANSFORM but because the
				// evolve hook is called from diagonalization which happens
				// well before the transform it doesn't
				// matter. Actually, it matters because using
				// BEFORE_TRANSFORM would return quantumNumbersOld
				// which aren't set before the changeOfBasis
				
				size_t qn = lrs_.super().qn(lrs_.super().partition(m));
				
				mettsStochastics_.update(qn,sites);
				
				size_t alphaFixed = mettsStochastics_.chooseRandomState(sites.first);
				size_t betaFixed = mettsStochastics_.chooseRandomState(sites.second);
				std::cerr<<"GETNEWPURES site="<<sites<<"\n";
				const MatrixType& transformSystem = 
				                         wft_.transform(ProgramGlobals::SYSTEM);
				VectorType newVector1(transformSystem.n_row());
				getNewPure(newVector1,pureVectors_.first,ProgramGlobals::SYSTEM,
				           alphaFixed,lrs_.left(),transformSystem,sites.first);
				pureVectors_.first = newVector1;
				
				const MatrixType& transformEnviron = 
				                        wft_.transform(ProgramGlobals::ENVIRON);
				VectorType newVector2(transformEnviron.n_row());
				getNewPure(newVector2,pureVectors_.second,ProgramGlobals::ENVIRON,
						   betaFixed,lrs_.right(),transformEnviron,sites.second);
				pureVectors_.second = newVector2;
				setFromInfinite(targetVectors_[0]);
			}

			void getFullVector(std::vector<RealType>& v,size_t m)
			{
				int offset = lrs_.super().partition(m);
				int total = lrs_.super().partition(m+1) - offset;

				PackIndicesType pack(lrs_.left().size());
				v.resize(total);
				for (int i=0;i<total;i++) {
					size_t alpha,beta;
					pack.unpack(alpha,beta,lrs_.super().permutation(i+offset));
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
				size_t ns = oldVector.size();
				size_t transformNrow =  (transform.n_row()==0) ? (ns*ns) : 
				                                           transform.n_row();
				newVector.resize(transformNrow);
				for (size_t gamma=0;gamma<transformNrow;gamma++) {
					newVector[gamma] = 0;
					for (size_t alpha=0;alpha<ns;alpha++) {
						size_t gammaPrime = (direction==ProgramGlobals::SYSTEM) ? 
						    basis.permutationInverse(alpha + alphaFixed*ns) :
						    basis.permutationInverse(alphaFixed + alpha*ns);
						if (transform.n_row()==0) {
							if (gamma == gammaPrime) 
								newVector[gamma] += oldVector[alpha];
							continue;
						}
						newVector[gamma] += transform(gamma,gammaPrime) * 
						                               oldVector[alpha];
					}
				}
			}
			
			void setInitialPure(VectorType& oldVector,size_t site)
			{
				size_t alphaFixed = mettsStochastics_.chooseRandomState(site);
				oldVector.resize(hilbertSizePerSite_);
				for (size_t i=0;i<oldVector.size();i++) {
					oldVector[i] = (i==alphaFixed) ? 1 : 0;
				}
			}
			
			void setFromInfinite(VectorWithOffsetType& phi)
			{
				phi = psi_;
				if (phi.sectors()!=1)
					throw std::runtime_error("MEttsTargetting::setFromInfinite(...): expected only one sector for g.s.\n");
				for (size_t ii=0;ii<phi.sectors();ii++) {
					size_t i0 = phi.sector(ii);
					VectorType v;
					getFullVector(v,i0);
					phi.setDataInSector(v,i0);
				}
			}

			void collapseVector(size_t m,
			                    VectorType& w,
			                    const VectorType& v,
			                    size_t alphaFixed,
			                    size_t betaFixed)
			{
				int offset = lrs_.super().partition(m);
				int total = lrs_.super().partition(m+1) - offset;

				size_t nk = hilbertSizePerSite_;
				size_t ns = lrs_.left().size();
				PackIndicesType packSuper(ns);
				PackIndicesType packLeft(ns/nk);
				PackIndicesType packRight(nk);
				for (size_t i=0;i<total;i++) {
					w[i] = 0;
					size_t alpha,beta;
					packSuper.unpack(alpha,beta,lrs_.super().permutation(i+offset));
					size_t alpha0,alpha1;
					packLeft.unpack(alpha0,alpha1,alpha);
					size_t beta0,beta1;
					packRight.unpack(beta0,beta1,beta);
					if (alpha1!=alphaFixed || beta0 != betaFixed) continue;
					w[i] = v[i];
				}
			}

			// in situ computation:
			void cocoon(size_t direction,const PairType& sites) const
			{
				std::cerr<<"-------------&*&*&* In-situ measurements start\n";
				test(psi_,psi_,direction,"<PSI|A|PSI>",sites);
				
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
				if (phi.size()==0) setFromInfinite(targetVectors_[startEnd.first]);
				
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
					// Only time differences here (i.e. betas_[i] not betas_[i]+currentBeta_)
					calcTargetVector(targetVectors_[i],
					        targetVectors_[startEnd.first],T,V,Eg,eigs,i,steps);
					//normalize(targetVectors_[i]);
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
				if (T.n_col()!=T.n_row())
					throw std::runtime_error("T is not square\n");
				if (V.n_col()!=T.n_col())
					throw std::runtime_error("V is not nxn2\n");
				// for (size_t j=0;j<v.size();j++) v[j] = 0; <-- harmful if v is sparse
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
				typename ModelType::ModelHelperType modelHelper(p,lrs_,model_.orbitals());
				 		//,useReflection_);
				typename LanczosSolverType::LanczosMatrixType lanczosHelper(&model_,&modelHelper);
			
				RealType eps= 0.01*ProgramGlobals::LanczosTolerance;
				size_t iter= ProgramGlobals::LanczosSteps;

				//srand48(3243447);
				LanczosSolverType lanczosSolver(lanczosHelper,iter,eps,parallelRank_);
				
				TridiagonalMatrixType ab;
				size_t total = phi.effectiveSize(i0);
				TargetVectorType phi2(total);
				phi.extract(phi2,i0);
				/* std::ostringstream msg;
				msg<<"Calling tridiagonalDecomposition...\n";
				progress_.printline(msg,std::cerr);*/
				lanczosSolver.tridiagonalDecomposition(phi2,ab,V);
				ab.buildDenseMatrix(T);
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
				if (V.n_col()>V.n_row()) throw std::runtime_error("cols > rows\n");
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

			void zeroOutVectors()
			{
				for (size_t i=0;i<targetVectors_.size();i++) 
					targetVectors_[i].resize(lrs_.super().size());
			}

//			void printHeader()
//			{
//				io_.print(mettsStruct_);
//				std::string label = "times";
//				io_.printVector(betas_,label);
//				label = "weights";
//				io_.printVector(weight_,label);
//				std::string s = "GsWeight="+ttos(gsWeight_);
//				io_.printline(s);
//			}

			void test(const VectorWithOffsetType& src1,
			          const VectorWithOffsetType& src2,
			          size_t systemOrEnviron,
			          const std::string& label,
			          const PairType& sites) const
			{
				throw std::runtime_error("Metts: test(...): not implemented\n");
// 				VectorWithOffsetType dest;
// 				OperatorType A = mettsStruct_.aOperators[0];
// 				PsimagLite::CrsMatrix<ComplexType> tmpC(model_.getOperator("c",0,0));
// 				PsimagLite::CrsMatrix<ComplexType> tmpCt;
// 				transposeConjugate(tmpCt,tmpC);
// 				multiply(A.data,tmpCt,tmpC);
// 				A.fermionSign = 1;
// 				//A.data = tmpC;
// 				FermionSign fs(lrs_.left(),mettsStruct_.electrons);
// 				applyOpLocal_(dest,src1,A,fs,systemOrEnviron);
// 
// 				ComplexType sum = 0;
// 				for (size_t ii=0;ii<dest.sectors();ii++) {
// 					size_t i = dest.sector(ii);
// 					size_t offset1 = dest.offset(i);
// 					for (size_t jj=0;jj<src2.sectors();jj++) {
// 						size_t j = src2.sector(jj);
// 						size_t offset2 = src2.offset(j);
// 						if (i!=j) continue; //throw std::runtime_error("Not same sector\n");
// 						for (size_t k=0;k<dest.effectiveSize(i);k++) 
// 							sum+= dest[k+offset1] * std::conj(src2[k+offset2]);
// 					}
// 				}
// 				std::cerr<<site<<" "<<sum<<" "<<" "<<currentBeta_;
// 				std::cerr<<" "<<label<<std::norm(src1)<<" "<<std::norm(src2)<<" "<<std::norm(dest)<<"\n";
			}


			size_t stage_;
			VectorWithOffsetType psi_;
			const LeftRightSuperType& lrs_;
			const ModelType& model_;
			const TargettingParamsType& mettsStruct_;
			const WaveFunctionTransfType& wft_;
			PsimagLite::ProgressIndicator progress_;
			RealType currentBeta_;
			std::vector<RealType> betas_,weight_;
			std::vector<VectorWithOffsetType> targetVectors_;
			RealType gsWeight_;
			//typename IoType::Out io_;
			ApplyOperatorType applyOpLocal_;
			size_t hilbertSizePerSite_;
			MettsStochasticsType mettsStochastics_;
			size_t timesWithoutAdvancement_;
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
