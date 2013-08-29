
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

#ifndef TIMESTEPTARGETTING_H
#define TIMESTEPTARGETTING_H

#include <iostream>
#include "ProgressIndicator.h"
#include "BLAS.h"
#include "ApplyOperatorLocal.h"
#include "TimeSerializer.h"
#include "TimeStepParams.h"
#include "ProgramGlobals.h"
#include "ParametersForSolver.h"
#include "ParallelWft.h"
#include "TimeVectorsKrylov.h"
#include "TimeVectorsRungeKutta.h"
#include "TimeVectorsSuzukiTrotter.h"
#include "CommonTargetting.h"

namespace Dmrg {
	template<template<typename,typename,typename> class LanczosSolverTemplate,
		 template<typename,typename> class InternalProductTemplate,
	         template<typename,typename> class WaveFunctionTransfTemplate,
	         typename ModelType_,
	         typename IoType_,
	         template<typename> class VectorWithOffsetTemplate>
	class TimeStepTargetting  {

		enum {BORDER_NEITHER, BORDER_LEFT, BORDER_RIGHT};

	public:

			typedef ModelType_ ModelType;
			typedef IoType_ IoType;
			typedef typename ModelType::RealType RealType;
			typedef std::complex<RealType> ComplexType;
			typedef InternalProductTemplate<ComplexType,ModelType>
				InternalProductType;
			typedef typename ModelType::OperatorsType OperatorsType;
			typedef typename ModelType::ModelHelperType ModelHelperType;
			typedef typename ModelHelperType::LeftRightSuperType
				LeftRightSuperType;
			typedef typename LeftRightSuperType::BasisWithOperatorsType
					BasisWithOperatorsType;
			typedef PsimagLite::Vector<SizeType>::Type VectorSizeType;
			//typedef BasisWithOperators<OperatorsType,ConcurrencyType> BasisWithOperatorsType;
			typedef typename PsimagLite::Vector<ComplexType>::Type ComplexVectorType;
			//typedef std::VectorWithOffset<ComplexType> VectorWithOffsetType;
			typedef PsimagLite::ParametersForSolver<RealType> ParametersForSolverType;
			typedef LanczosSolverTemplate<ParametersForSolverType,InternalProductType,ComplexVectorType> LanczosSolverType;
			typedef typename PsimagLite::Vector<RealType>::Type VectorType;
			//typedef typename BasisWithOperatorsType::SparseMatrixType SparseMatrixType;
			typedef PsimagLite::Matrix<ComplexType> ComplexMatrixType;
			typedef typename BasisWithOperatorsType::OperatorType OperatorType;
			typedef typename BasisWithOperatorsType::BasisType BasisType;
			typedef TimeStepParams<ModelType> TargettingParamsType;
			typedef typename BasisType::BlockType BlockType;
			typedef VectorWithOffsetTemplate<ComplexType> VectorWithOffsetType;
			typedef WaveFunctionTransfTemplate<LeftRightSuperType,VectorWithOffsetType> WaveFunctionTransfType;
			typedef ComplexVectorType TargetVectorType;
			typedef BlockMatrix<ComplexMatrixType> ComplexBlockMatrixType;
			typedef ApplyOperatorLocal<LeftRightSuperType,VectorWithOffsetType,TargetVectorType> ApplyOperatorType;
			typedef TimeSerializer<VectorWithOffsetType> TimeSerializerType;
			typedef typename OperatorType::SparseMatrixType SparseMatrixType;
			typedef typename BasisWithOperatorsType::BasisDataType BasisDataType;
			typedef TimeVectorsBase<TargettingParamsType,ModelType,WaveFunctionTransfType,
									LanczosSolverType,VectorWithOffsetType> TimeVectorsBaseType;
			typedef TimeVectorsKrylov<TargettingParamsType,ModelType,WaveFunctionTransfType,
									  LanczosSolverType,VectorWithOffsetType> TimeVectorsKrylovType;
			typedef TimeVectorsRungeKutta<TargettingParamsType,ModelType,WaveFunctionTransfType,
										  LanczosSolverType,VectorWithOffsetType> TimeVectorsRungeKuttaType;
			typedef TimeVectorsSuzukiTrotter<TargettingParamsType,ModelType,WaveFunctionTransfType,
											 LanczosSolverType,VectorWithOffsetType> TimeVectorsSuzukiTrotterType;
			typedef CommonTargetting<ModelType,TargettingParamsType,WaveFunctionTransfType,VectorWithOffsetType,LanczosSolverType>
		        CommonTargettingType;
			typedef typename CommonTargettingType::VectorStringType VectorStringType;

			enum {DISABLED,OPERATOR,WFT_NOADVANCE,WFT_ADVANCE};
			enum {EXPAND_ENVIRON=WaveFunctionTransfType::EXPAND_ENVIRON,
			EXPAND_SYSTEM=WaveFunctionTransfType::EXPAND_SYSTEM,
			INFINITE=WaveFunctionTransfType::INFINITE};

			static SizeType const PRODUCT = TargettingParamsType::PRODUCT;
			static SizeType const SUM = TargettingParamsType::SUM;

			TimeStepTargetting(const LeftRightSuperType& lrs,
	 		                   const ModelType& model,
			                   const TargettingParamsType& tstStruct,
			                   const WaveFunctionTransfType& wft,
			                   const SizeType& quantumSector) // quantumSector is ignored here
			: stage_(tstStruct.sites.size(),DISABLED),
			  lrs_(lrs),
			  model_(model),
			  tstStruct_(tstStruct),
			  wft_(wft),
			  progress_("TimeStepTargetting"),
			  currentTime_(0),
			  times_(tstStruct_.timeSteps),
			  weight_(tstStruct_.timeSteps),
			  targetVectors_(tstStruct_.timeSteps),
			  commonTargetting_(lrs,model,tstStruct),
			  applyOpLocal_(lrs),
			  E0_(0),
			  timeVectorsBase_(0)
			{
				if (!wft.isEnabled()) throw PsimagLite::RuntimeError
				       (" TimeStepTargetting needs an enabled wft\n");
				if (tstStruct_.sites.size()==0) throw PsimagLite::RuntimeError
				       (" TimeStepTargetting needs at least one TSPSite\n");

				RealType tau =tstStruct_.tau;
				RealType sum = 0;
				SizeType n = times_.size();
				gsWeight_ = (tstStruct_.concatenation == SUM) ? 0.1 : 0.0;
				RealType factor = (n+4.0)/(n+2.0);
				factor *= (1.0 - gsWeight_);
				for (SizeType i=0;i<n;i++) {
					times_[i] = i*tau/(n-1);
					weight_[i] = factor/(n+4);
					sum += weight_[i];
				}
				sum -= weight_[0];
				sum -= weight_[n-1];
				weight_[0] = weight_[n-1] = 2*factor/(n+4);
				sum += weight_[n-1];
				sum += weight_[0];
				
				gsWeight_=1.0-sum;
				sum += gsWeight_;
				assert(fabs(sum-1.0)<1e-5);
//				std::cerr<<"GSWEIGHT="<<gsWeight_<<"\n";

				PsimagLite::String s (__FILE__);
				s += " Unknown algorithm\n";
				switch (tstStruct_.algorithm) {
				case TargettingParamsType::KRYLOV:
					timeVectorsBase_ = new TimeVectorsKrylovType(
								currentTime_,tstStruct_,times_,targetVectors_,model_,wft_,lrs_,E0_);
					break;
				case TargettingParamsType::RUNGE_KUTTA:
					timeVectorsBase_ = new TimeVectorsRungeKuttaType(
								currentTime_,tstStruct_,times_,targetVectors_,model_,wft_,lrs_,E0_);
					break;
				case TargettingParamsType::SUZUKI_TROTTER:
					timeVectorsBase_ = new TimeVectorsSuzukiTrotterType(
								currentTime_,tstStruct_,times_,targetVectors_,model_,wft_,lrs_,E0_,&nonZeroQns_);
					break;
				default:
					throw PsimagLite::RuntimeError(s.c_str());
				}
			}

			~TimeStepTargetting()
			{
				if (timeVectorsBase_)
					delete timeVectorsBase_;
			}

			const ModelType& model() const { return model_; }

			RealType weight(SizeType i) const
			{
				assert(!commonTargetting_.allStages(DISABLED,stage_));
				return weight_[i];
			}

			RealType gsWeight() const
			{
				if (commonTargetting_.allStages(DISABLED,stage_)) return 1.0;
				return gsWeight_;
			}

			RealType normSquared(SizeType i) const
			{
				return commonTargetting_.normSquared(targetVectors_[i]);
			}

			template<typename SomeBasisType>
			void setGs(const typename PsimagLite::Vector<TargetVectorType>::Type& v,
				   const SomeBasisType& someBasis)
			{
				psi_.set(v,someBasis);
				assert(psi_.size()==lrs_.super().size());
			}

			const ComplexType& operator[](SizeType i) const { return psi_[i]; }

			ComplexType& operator[](SizeType i) { return psi_[i]; }

			const VectorWithOffsetType& gs() const { return psi_; }

			bool includeGroundStage() const
			{
				if (!commonTargetting_.noStageIs(DISABLED,stage_)) return true;
				bool b = (fabs(gsWeight_)>1e-6);

//				std::cerr<<"includeGroundState="<<b<<"\n";
				return b;
			}

			const RealType& time() const {return currentTime_; }

			SizeType size() const
			{
				if (commonTargetting_.allStages(DISABLED,stage_)) return 0;
				return targetVectors_.size();
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

				VectorWithOffsetType phiOld = psi_;
				VectorWithOffsetType phiNew;
				VectorWithOffsetType vectorSum;
				SizeType max = tstStruct_.sites.size();

				if (commonTargetting_.noStageIs(DISABLED,stage_)) {
					max = 1;
					for (SizeType i=0;i<stage_.size();i++) {
						if (stage_[i]==OPERATOR) stage_[i] = WFT_NOADVANCE;
						if (stage_[i]==WFT_ADVANCE) stage_[i] = WFT_NOADVANCE;
					}
				}

				// Loop over each operator that needs to be applied 
				// in turn to the g.s.

				for (SizeType i=0;i<max;i++) {
					bool wasAnOperatorApplied = false;
					if (!evolve(i,phiNew,phiOld,Eg,direction,block1,loopNumber,max-1,wasAnOperatorApplied)) {
						continue;
					}
					if (tstStruct_.concatenation==PRODUCT) {
						phiOld = phiNew;
					} else {
						if (wasAnOperatorApplied) vectorSum += phiNew;
						else vectorSum = phiNew;
					}
				}

				//std::cerr<<"site="<<block1[0]<<" COUNT="<<count<<"\n";
				if (tstStruct_.concatenation==SUM) phiNew = vectorSum;

				typename TimeVectorsBaseType::PairType startEnd(0,times_.size());
				bool allOperatorsApplied = (commonTargetting_.noStageIs(DISABLED,stage_) &&
				                            commonTargetting_.noStageIs(OPERATOR,stage_));

				timeVectorsBase_->calcTimeVectors(startEnd,
				                                  Eg,
				                                  phiNew,
				                                  direction,
				                                  allOperatorsApplied,
				                                  block1);
				
				cocoon(direction,block1); // in-situ
				printEnergies(); // in-situ
				printNormsAndWeights();
			}

			void load(const PsimagLite::String& f)
			{
				for (SizeType i=0;i<stage_.size();i++) stage_[i] = WFT_NOADVANCE;

				typename IoType::In io(f);

				TimeSerializerType ts(io,IoType::In::LAST_INSTANCE);
				for (SizeType i=0;i<targetVectors_.size();i++) targetVectors_[i] = ts.vector(i);
				currentTime_ = ts.time();

				psi_.load(io,"PSI");
			}

			void print(std::ostream& os) const
			{
				os<<"TSTWeightsTimeVectors=";
				for (SizeType i=0;i<weight_.size();i++)
					os<<weight_[i]<<" ";
				os<<"\n";
				os<<"TSTWeightGroundState="<<gsWeight_<<"\n";
			}

			bool evolve(SizeType i,
				    VectorWithOffsetType& phiNew,
				    const VectorWithOffsetType& phiOld,
				    RealType Eg,
				    SizeType direction,
				    const BlockType& block,
				    SizeType loopNumber,
				    SizeType lastI,
				    bool& wasAnOperatorApplied)
			{
				static SizeType  timesWithoutAdvancement=0;
				static bool firstSeeLeftCorner = false;
				wasAnOperatorApplied = false;

				if (direction==INFINITE) {
					E0_ = Eg;
					return false;
				}
				if (tstStruct_.startingLoops.size()>0 && tstStruct_.startingLoops[i]>loopNumber)
					return false;

				assert(block.size()==1);
				SizeType site = block[0];

				if (site != tstStruct_.sites[i] && stage_[i]==DISABLED)
					return false;

				if (site != tstStruct_.sites[i] && stage_[i]!=DISABLED && i>0)
					return false;

				if (site == tstStruct_.sites[i] && stage_[i]==DISABLED) {
					stage_[i]=OPERATOR;
					wasAnOperatorApplied=true;
				} else {
					stage_[i]=WFT_NOADVANCE;
				}
				if (stage_[i] == OPERATOR) commonTargetting_.checkOrder(i,stage_);

				if (timesWithoutAdvancement >= tstStruct_.advanceEach) {
					stage_[i] = WFT_ADVANCE;
					if (i==lastI) {
						currentTime_ += tstStruct_.tau;
						timesWithoutAdvancement=1;
					}
				} else {
					if (i==lastI && stage_[i]==WFT_NOADVANCE && firstSeeLeftCorner)
						timesWithoutAdvancement++;
				}

				if (!firstSeeLeftCorner && i==lastI && stage_[i]==WFT_NOADVANCE && site==1)
					firstSeeLeftCorner=true;

				PsimagLite::OstringStream msg2;
				msg2<<"Steps without advance: "<<timesWithoutAdvancement<<" site="<<site<<" currenTime="<<currentTime_;
				if (timesWithoutAdvancement>0) progress_.printline(msg2,std::cout);
				
				PsimagLite::OstringStream msg;
				msg<<"Evolving, stage="<<commonTargetting_.getStage(i,stage_);
				msg<<" site="<<site<<" loopNumber="<<loopNumber;
				msg<<" Eg="<<Eg;
				progress_.printline(msg,std::cout);
				
				// phi = A|psi>
				computePhi(i,phiNew,phiOld,direction,block);
				
				return true;
			}

			void initialGuess(VectorWithOffsetType& v,const VectorSizeType& block) const
			{
				VectorSizeType nk;
				commonTargetting_.setNk(nk,block);
				wft_.setInitialVector(v,psi_,lrs_,nk);
			}

			template<typename IoOutputType>
			void save(const VectorSizeType& block,IoOutputType& io) const
			{
				PsimagLite::OstringStream msg;
				msg<<"Saving state...";
				progress_.printline(msg,std::cout);

				SizeType marker = 0;
				if (commonTargetting_.noStageIs(DISABLED,stage_)) marker = 1;

				TimeSerializerType ts(currentTime_,block[0],targetVectors_,marker);
				ts.save(io);
				psi_.save(io,"PSI");
			}

			void updateOnSiteForTimeDep(BasisWithOperatorsType& basisWithOps) const
			{

				BlockType X = basisWithOps.block();
				if (X.size()!=1) return;
				assert(X[0]==0 || X[0]==lrs_.super().block().size()-1);
				typename PsimagLite::Vector<OperatorType>::Type creationMatrix;
				SparseMatrixType hmatrix;
				BasisDataType q;
				model_.setNaturalBasis(creationMatrix,hmatrix,q,X,currentTime_);
				basisWithOps.setVarious(X,hmatrix,q,creationMatrix);
			}

			bool end() const
			{
				return (tstStruct_.maxTime != 0 && currentTime_ >= tstStruct_.maxTime);
			}

		private:

			void printNormsAndWeights() const
			{
				if (commonTargetting_.allStages(DISABLED,stage_)) return;

				PsimagLite::OstringStream msg;
				msg<<"gsWeight="<<gsWeight_<<" weights= ";
				for (SizeType i = 0; i < weight_.size(); i++)
					msg<<weight_[i]<<" ";
				progress_.printline(msg,std::cout);

				PsimagLite::OstringStream msg2;
				msg2<<"gsNorm="<<std::norm(psi_)<<" norms= ";
				for (SizeType i = 0; i < weight_.size(); i++)
					msg2<<normSquared(i)<<" ";
				progress_.printline(msg2,std::cout);
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
				msg<<"Hamiltonian average at time="<<currentTime_<<" for target="<<whatTarget;
				msg<<" sector="<<i0<<" <phi(t)|H|phi(t)>="<<(phi2*x)<<" <phi(t)|phi(t)>="<<(phi2*phi2);
				progress_.printline(msg,std::cout);
			}

			// in situ computation:
			void cocoon(SizeType direction,const BlockType& block) const
			{
				std::cout<<"-------------&*&*&* In-situ measurements start\n";

				if (commonTargetting_.noStageIs(DISABLED,stage_))
					std::cout<<"ALL OPERATORS HAVE BEEN APPLIED\n";
				else
					std::cout<<"NOT ALL OPERATORS APPLIED YET\n";

				PsimagLite::String modelName = model_.params().model;

				if (modelName == "HubbardOneBand" ||
				    modelName == "HubbardOneBandExtended") {
					cocoon_(direction,block);
				}

				SizeType site = block[0];

				commonTargetting_.cocoon(direction,site,psi_,"PSI");
				commonTargetting_.cocoon(direction,site,targetVectors_[0],"P0");

				std::cout<<"-------------&*&*&* In-situ measurements end\n";
			}


			void cocoon_(SizeType direction,const BlockType& block) const
			{
				SizeType site = block[0];
				PsimagLite::CrsMatrix<ComplexType> tmpC(model_.naturalOperator("nup",0,0));
				//PsimagLite::CrsMatrix<ComplexType> tmpCt;
				//transposeConjugate(tmpCt,tmpC);
				//multiply(A.data,tmpCt,tmpC);
				int fermionSign1 = 1;
				const std::pair<SizeType,SizeType> jm1(0,0);
				RealType angularFactor1 = 1.0;
				typename OperatorType::Su2RelatedType su2Related1;
				OperatorType nup(tmpC,fermionSign1,jm1,angularFactor1,su2Related1);

				nup.data = tmpC;
				nup.fermionSign = 1;
				//A.data = tmpC;

				test(psi_,psi_,direction,"<PSI|nup|PSI>",site,nup);
				PsimagLite::String s = "<P0|nup|P0>";
				test(targetVectors_[0],targetVectors_[0],direction,s,site,nup);

				PsimagLite::CrsMatrix<ComplexType> tmpC2(model_.naturalOperator("ndown",0,0));
				OperatorType ndown(tmpC2,fermionSign1,jm1,angularFactor1,su2Related1);
				//if (includeGroundStage())
				test(psi_,psi_,direction,"<PSI|ndown|PSI>",site,ndown);
				s = "<P0|ndown|P0>";
				test(targetVectors_[0],targetVectors_[0],direction,s,site,ndown);

				PsimagLite::CrsMatrix<ComplexType> tmpC3 = (nup.data * ndown.data);
				OperatorType doubleOcc(tmpC3,fermionSign1,jm1,angularFactor1,su2Related1);
				test(psi_,psi_,direction,"<PSI|doubleOcc|PSI>",site,doubleOcc);
				s = "<P0|doubleOcc|P0>";
				test(targetVectors_[0],targetVectors_[0],direction,s,site,doubleOcc);
			}

			void computePhi(SizeType i,
			                VectorWithOffsetType& phiNew,
			                const VectorWithOffsetType& phiOld,
			                SizeType systemOrEnviron,
			                const VectorSizeType& block)
			{
				if (block.size()!=1) {
					PsimagLite::String str(__FILE__);
					str += " " + ttos(__LINE__) + "\n";
					str += "computePhi only blocks of one site supported\n";
					throw PsimagLite::RuntimeError(str.c_str());
				}
				VectorSizeType nk;
				commonTargetting_.setNk(nk,block);
				SizeType site = block[0];

				SizeType indexAdvance = times_.size()-1;
				SizeType indexNoAdvance = 0;
				if (stage_[i]==OPERATOR) {
					PsimagLite::OstringStream msg;
					msg<<"I'm applying a local operator at site "<<site<<" now";
					progress_.printline(msg,std::cout);
					VectorSizeType electrons;
					model_.findElectronsOfOneSite(electrons,site);
					FermionSign fs(lrs_.left(),electrons);
					applyOpLocal_(phiNew,phiOld,tstStruct_.aOperators[i],fs,systemOrEnviron);

					applyOperatorAtBorder(phiNew,phiOld,systemOrEnviron,block);

					setQuantumNumbers(phiNew);

				} else if (stage_[i]== WFT_NOADVANCE || stage_[i]== WFT_ADVANCE) {
					SizeType advance = indexNoAdvance;
					if (stage_[i] == WFT_ADVANCE) {
						advance = indexAdvance;
						timeVectorsBase_->timeHasAdvanced();
					}
					PsimagLite::OstringStream msg;
					msg<<"I'm calling the WFT now";
					progress_.printline(msg,std::cout);

					if (tstStruct_.aOperators.size()==1)
						guessPhiSectors(phiNew,i,systemOrEnviron,site);
					else
						phiNew.populateFromQns(nonZeroQns_,lrs_.super());

					// OK, now that we got the partition number right, let's wft:
					wft_.setInitialVector(phiNew,targetVectors_[advance],lrs_,nk); // generalize for su(2)
//					phiNew.collapseSectors();
//					std::cerr<<"WFT --> NORM of phiNew="<<norm(phiNew)<<" NORM of tv="<<norm(targetVectors_[advance])<<" when i="<<i<<" advance="<<advance<<"\n";

				} else {
					throw PsimagLite::RuntimeError("It's 5 am, do you know what line "
						" your code is exec-ing?\n");
				}
				RealType norma = norm(phiNew);
				if (norma==0) throw PsimagLite::RuntimeError("Norm of phi is zero\n");

			}

			void applyOperatorAtBorder(VectorWithOffsetType& phiNew,
			                           const VectorWithOffsetType& phiOld,
			                           SizeType systemOrEnviron,
			                           const VectorSizeType& block)
			{
				SizeType whichBorder = findWhichBorder(block);
				if (whichBorder == BORDER_NEITHER) return;

				int iOfBorder = findBorderIndex(whichBorder);
				if (iOfBorder < 0) return;

				if (block.size() != 1) {
					PsimagLite::String str("applyOperatorAtBorder(): ");
					str += " implemented only for block.size() == 1\n";
					throw PsimagLite::RuntimeError(str);
				}

				SizeType site = tstStruct_.sites[iOfBorder];

				stage_[iOfBorder] = OPERATOR;

				PsimagLite::OstringStream msg;
				msg<<"I'm applying a local operator at site "<<site<<" now";
				progress_.printline(msg,std::cout);
				VectorSizeType electrons;
				model_.findElectronsOfOneSite(electrons,site);
				FermionSign fs(lrs_.left(),electrons);
				VectorWithOffsetType phiNew2;
				const VectorWithOffsetType& phiSrc = (tstStruct_.concatenation==PRODUCT) ?
				            phiNew : phiOld;

				applyOpLocal_(phiNew2,
				              phiSrc,
				              tstStruct_.aOperators[iOfBorder],
				              fs,
				              systemOrEnviron,
				              true);

				if (tstStruct_.concatenation==PRODUCT)
					phiNew = phiNew2;
				else
					phiNew += phiNew2;
			}

			SizeType findWhichBorder(const VectorSizeType& block) const
			{
				SizeType sitesPerBlock = model_.params().sitesPerBlock;

				SizeType siteMin = *std::min_element(block.begin(),block.end());
				if (siteMin <= sitesPerBlock) return BORDER_LEFT;

				SizeType siteMax = *std::max_element(block.begin(),block.end());
				SizeType linSize = model_.geometry().numberOfSites();

				if (siteMax + sitesPerBlock + 1 >= linSize) return BORDER_RIGHT;

				return BORDER_NEITHER;

			}

			int findBorderIndex(SizeType whichBorder) const
			{
				assert(whichBorder == BORDER_LEFT || whichBorder == BORDER_RIGHT);

				SizeType sitesPerBlock = model_.params().sitesPerBlock;
				SizeType linSize = model_.geometry().numberOfSites();
				assert(sitesPerBlock < linSize);
				SizeType start = linSize - sitesPerBlock;

				for (SizeType j = 0; j < tstStruct_.sites.size(); ++j) {
					if (whichBorder == BORDER_LEFT) {
						for (SizeType i = 0; i < sitesPerBlock; ++i) {
							if (tstStruct_.sites[j] == i) return j;
						}
					} else {
						for (SizeType i = start; i < linSize; ++i) {
							if (tstStruct_.sites[j] == i) return j;
						}
					}
				}

				return -1;
			}

			void checkNorms() const
			{
				PsimagLite::OstringStream msg;
				msg<<"Checking norms: ";
				for (SizeType i=0;i<targetVectors_.size();i++) {
					RealType norma = std::norm(targetVectors_[i]);
					msg<<" norma["<<i<<"]="<<norma;
					assert(norma>1e-10);
				}
				progress_.printline(msg,std::cout);
			}



			void check1(const VectorWithOffsetType& phi,SizeType i0) const
			{
				ComplexType ret = 0;
				SizeType total = phi.effectiveSize(i0);

				for (SizeType j=0;j<total;j++)
					ret += std::conj(phi.fastAccess(i0,j))*phi.fastAccess(i0,j);
				std::cerr<<"check1, norm  of phi="<<ret<<"\n";

			}

			void check2(const ComplexMatrixType& T,
				    const ComplexMatrixType& V,
				    const VectorWithOffsetType& phi,
				    SizeType n2,
				    SizeType i0) const
			{
				check3(V);

				ComplexVectorType r(n2);

				for (SizeType k=0;k<n2;k++) {
					ComplexType sum = 0.0;
					for (SizeType kprime=0;kprime<n2;kprime++) {
						if (kprime!=k) continue;
						ComplexType tmpV = calcVTimesPhi(kprime,V,phi,i0);
//						sum += conj(T(kprime,k))*tmpV;
						sum += tmpV;
					}
					r[k] = sum;
				}
				std::cerr<<"check2, norm  of  V . phi="<<PsimagLite::norm(r)<<"\n";

			}

			void check3(const ComplexMatrixType& V) const
			{
				for (SizeType i=0;i<V.n_col();i++) {
					ComplexType sum = 0;
					for (SizeType j=0;j<V.n_row();j++) {
						sum += std::conj(V(j,i)) * V(j,0);
					}
					std::cerr<<"V["<<i<<"] * V[0] = "<<sum<<"\n";
				}
			}

			//! This check is invalid if there are more than one sector
			void check1(const ComplexMatrixType& V,const TargetVectorType& phi2)
			{
				assert(V.n_col()<=V.n_row());
				TargetVectorType r(V.n_col());
				for (SizeType k=0;k<V.n_col();k++) {
					r[k] = 0.0;
					for (SizeType j=0;j<V.n_row();j++) 
						r[k] += conj(V(j,k))*phi2[j];
					// is r(k) == \delta(k,0)
					if (k==0 && std::norm(r[k]-1.0)>1e-5) 
						std::cerr<<"WARNING: r[0]="<<r[0]<<" != 1\n";
					if (k>0 && std::norm(r[k])>1e-5) 
						std::cerr<<"WARNING: r["<<k<<"]="<<r[k]<<" !=0\n";
				}
			}

			void setQuantumNumbers(const VectorWithOffsetType& v)
			{
				nonZeroQns_.clear();
				assert(v.size()==lrs_.super().size());
				for (SizeType i=0;i<v.sectors();i++) {
					SizeType i0 = v.sector(i);
					SizeType state = v.offset(i0);
					SizeType qn = lrs_.super().qn(state);
					nonZeroQns_.push_back(qn);
				}
			}

			void guessPhiSectors(VectorWithOffsetType& phi,SizeType i,SizeType systemOrEnviron,SizeType site)
			{
				VectorSizeType electrons;
				model_.findElectronsOfOneSite(electrons,site);
				FermionSign fs(lrs_.left(),electrons);
				if (commonTargetting_.allStages(WFT_NOADVANCE,stage_)) {
					VectorWithOffsetType tmpVector = psi_;
					for (SizeType j=0;j<tstStruct_.aOperators.size();j++) {
						applyOpLocal_(phi,tmpVector,tstStruct_.aOperators[j],fs,
							systemOrEnviron);
						tmpVector = phi;
					}
					return;
				}
				applyOpLocal_(phi,psi_,tstStruct_.aOperators[i],fs,systemOrEnviron);
			}

			void test(const VectorWithOffsetType& src1,
				  const VectorWithOffsetType& src2,
				  SizeType systemOrEnviron,
				  const PsimagLite::String& label,
				  SizeType site,
				  const OperatorType& A) const
			{
				VectorSizeType electrons;
				SizeType lastIndex = lrs_.left().block().size();
				assert(lastIndex>0);
				lastIndex--;
				SizeType siteCorrected = lrs_.left().block()[lastIndex];
				model_.findElectronsOfOneSite(electrons,siteCorrected);
				FermionSign fs(lrs_.left(),electrons);
				VectorWithOffsetType dest;
				applyOpLocal_(dest,src1,A,fs,systemOrEnviron);

				ComplexType sum = 0;
				for (SizeType ii=0;ii<dest.sectors();ii++) {
					SizeType i = dest.sector(ii);
					SizeType offset1 = dest.offset(i);
					for (SizeType jj=0;jj<src2.sectors();jj++) {
						SizeType j = src2.sector(jj);
						SizeType offset2 = src2.offset(j);
						if (i!=j) continue; //throw PsimagLite::RuntimeError("Not same sector\n");
						for (SizeType k=0;k<dest.effectiveSize(i);k++) 
							sum+= dest[k+offset1] * conj(src2[k+offset2]);
					}
				}
				std::cout<<site<<" "<<sum<<" "<<" "<<currentTime_;
				std::cout<<" "<<label<<" "<<(src1*src2)<<"\n";
			}

			VectorSizeType stage_;
			VectorWithOffsetType psi_;
			const LeftRightSuperType& lrs_;
			const ModelType& model_;
			const TargettingParamsType& tstStruct_;
			const WaveFunctionTransfType& wft_;
			PsimagLite::ProgressIndicator progress_;
			RealType currentTime_;
			typename PsimagLite::Vector<RealType>::Type times_,weight_;
			typename PsimagLite::Vector<VectorWithOffsetType>::Type targetVectors_;
			CommonTargettingType commonTargetting_;
			RealType gsWeight_;
			//typename IoType::Out io_;
			ApplyOperatorType applyOpLocal_;
			RealType E0_;
			TimeVectorsBaseType* timeVectorsBase_;
			VectorSizeType nonZeroQns_;

	};     //class TimeStepTargetting

	template<
		template<typename,typename,typename> class LanczosSolverTemplate,
			template<typename,typename> class InternalProductTemplate,
 		template<typename,typename> class WaveFunctionTransfTemplate,
			typename ModelType_,
			typename IoType_,
   			template<typename> class VectorWithOffsetTemplate>
	std::ostream& operator<<(std::ostream& os,
			const TimeStepTargetting<LanczosSolverTemplate,
			InternalProductTemplate,
			WaveFunctionTransfTemplate,ModelType_,IoType_,
			VectorWithOffsetTemplate>& tst)
	{
		tst.print(os);
		return os;
	}
} // namespace Dmrg

#endif
