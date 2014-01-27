
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
         typename MatrixVectorType_,
         typename WaveFunctionTransfType_,
         typename IoType_>
	class TimeStepTargetting  {

		enum {BORDER_NEITHER, BORDER_LEFT, BORDER_RIGHT};

	public:

		typedef MatrixVectorType_ MatrixVectorType;
			typedef typename MatrixVectorType::ModelType ModelType;
			typedef IoType_ IoType;
			typedef typename ModelType::RealType RealType;
			typedef std::complex<RealType> ComplexType;
			typedef typename ModelType::OperatorsType OperatorsType;
			typedef typename ModelType::ModelHelperType ModelHelperType;
			typedef typename ModelHelperType::LeftRightSuperType
				LeftRightSuperType;
			typedef typename LeftRightSuperType::BasisWithOperatorsType
					BasisWithOperatorsType;
			typedef PsimagLite::Vector<SizeType>::Type VectorSizeType;
			typedef WaveFunctionTransfType_ WaveFunctionTransfType;
			typedef typename WaveFunctionTransfType::VectorWithOffsetType VectorWithOffsetType;
			typedef typename VectorWithOffsetType::VectorType TargetVectorType;
			typedef PsimagLite::ParametersForSolver<RealType> ParametersForSolverType;
			typedef LanczosSolverTemplate<ParametersForSolverType,MatrixVectorType,TargetVectorType> LanczosSolverType;
			typedef typename PsimagLite::Vector<RealType>::Type VectorType;
			typedef PsimagLite::Matrix<ComplexType> ComplexMatrixType;
			typedef typename BasisWithOperatorsType::OperatorType OperatorType;
			typedef typename BasisWithOperatorsType::BasisType BasisType;
			typedef TimeStepParams<ModelType> TargettingParamsType;
			typedef typename BasisType::BlockType BlockType;
			typedef BlockMatrix<ComplexMatrixType> ComplexBlockMatrixType;
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
		typedef TargetHelper<ModelType,
		                     TargettingParamsType,
		                     WaveFunctionTransfType> TargetHelperType;
		typedef CommonTargetting<TargetHelperType,
		                         VectorWithOffsetType,
		                         LanczosSolverType> CommonTargettingType;
		typedef typename CommonTargettingType::VectorStringType VectorStringType;
		typedef typename CommonTargettingType::VectorVectorWithOffsetType VectorVectorWithOffsetType;

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
			: lrs_(lrs),
			  model_(model),
			  tstStruct_(tstStruct),
			  wft_(wft),
			  progress_("TimeStepTargetting"),
			  times_(tstStruct_.timeSteps),
			  weight_(tstStruct_.timeSteps),
			  commonTargetting_(lrs,model,tstStruct,wft)
			{
				if (!wft.isEnabled()) throw PsimagLite::RuntimeError
				       (" TimeStepTargetting needs an enabled wft\n");
				if (tstStruct_.sites.size()==0) throw PsimagLite::RuntimeError
				       (" TimeStepTargetting needs at least one TSPSite\n");

				RealType tau =tstStruct_.tau;
				RealType sum = 0;
				SizeType n = times_.size();
				gsWeight_ = (tstStruct_.concatenation == SUM) ? 0.1 : 0.0;
				gsWeight_ = commonTargetting_.setGsWeight(gsWeight_);

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

				commonTargetting_.initTimeVectors(times_);
			}

			const ModelType& model() const { return model_; }

			RealType weight(SizeType i) const
			{
				assert(!commonTargetting_.allStages(DISABLED));
				return weight_[i];
			}

			RealType gsWeight() const
			{
				if (commonTargetting_.allStages(DISABLED)) return 1.0;
				return gsWeight_;
			}

			RealType normSquared(SizeType i) const
			{
				return commonTargetting_.normSquared(i);
			}

			template<typename SomeBasisType>
			void setGs(const typename PsimagLite::Vector<TargetVectorType>::Type& v,
				   const SomeBasisType& someBasis)
			{
				commonTargetting_.psi().set(v,someBasis);
				assert(commonTargetting_.psi().size()==lrs_.super().size());
			}

			const ComplexType& operator[](SizeType i) const { return commonTargetting_.psi()[i]; }

			ComplexType& operator[](SizeType i) { return commonTargetting_.psi()[i]; }

			const VectorWithOffsetType& gs() const { return commonTargetting_.psi(); }

			bool includeGroundStage() const
			{
				if (!commonTargetting_.noStageIs(DISABLED)) return true;
				bool b = (fabs(gsWeight_)>1e-6);

//				std::cerr<<"includeGroundState="<<b<<"\n";
				return b;
			}

			const RealType& time() const {return commonTargetting_.currentTime(); }

			SizeType size() const
			{
				if (commonTargetting_.allStages(DISABLED)) return 0;
				return commonTargetting_.targetVectors().size();
			}

			const VectorWithOffsetType& operator()(SizeType i) const
			{
				return commonTargetting_.targetVectors()[i];
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
				VectorWithOffsetType phiNew;
				commonTargetting_.getPhi(phiNew,Eg,direction,block1[0],loopNumber);

				typename TimeVectorsBaseType::PairType startEnd(0,times_.size());
				bool allOperatorsApplied = (commonTargetting_.noStageIs(DISABLED) &&
				                            commonTargetting_.noStageIs(OPERATOR));

				commonTargetting_.calcTimeVectors(startEnd,
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
				commonTargetting_.setAllStagesTo(WFT_NOADVANCE);

				typename IoType::In io(f);

				TimeSerializerType ts(io,IoType::In::LAST_INSTANCE);
				for (SizeType i=0;i<commonTargetting_.targetVectors().size();i++)
					commonTargetting_.targetVectors(i) = ts.vector(i);

				commonTargetting_.setTime(ts.time());

				commonTargetting_.psi().load(io,"PSI");
			}

			void print(std::ostream& os) const
			{
				os<<"TSTWeightsTimeVectors=";
				for (SizeType i=0;i<weight_.size();i++)
					os<<weight_[i]<<" ";
				os<<"\n";
				os<<"TSTWeightGroundState="<<gsWeight_<<"\n";
			}

			void initialGuess(VectorWithOffsetType& v,const VectorSizeType& block) const
			{
				commonTargetting_.initialGuess(v,block);
			}

			template<typename IoOutputType>
			void save(const VectorSizeType& block,IoOutputType& io) const
			{
				PsimagLite::OstringStream msg;
				msg<<"Saving state...";
				progress_.printline(msg,std::cout);

				SizeType marker = 0;
				if (commonTargetting_.noStageIs(DISABLED)) marker = 1;

				TimeSerializerType ts(commonTargetting_.currentTime(),
				                      block[0],
				                      commonTargetting_.targetVectors(),
				                      marker);
				ts.save(io);
				commonTargetting_.psi().save(io,"PSI");
			}

			void updateOnSiteForTimeDep(BasisWithOperatorsType& basisWithOps) const
			{

				BlockType X = basisWithOps.block();
				if (X.size()!=1) return;
				assert(X[0]==0 || X[0]==lrs_.super().block().size()-1);
				typename PsimagLite::Vector<OperatorType>::Type creationMatrix;
				SparseMatrixType hmatrix;
				BasisDataType q;
				model_.setNaturalBasis(creationMatrix,hmatrix,q,X,commonTargetting_.currentTime());
				basisWithOps.setVarious(X,hmatrix,q,creationMatrix);
			}

			bool end() const
			{
				return (tstStruct_.maxTime != 0 &&
				        commonTargetting_.currentTime() >= tstStruct_.maxTime);
			}

		private:

			void printNormsAndWeights() const
			{
				if (commonTargetting_.allStages(DISABLED)) return;

				PsimagLite::OstringStream msg;
				msg<<"gsWeight="<<gsWeight_<<" weights= ";
				for (SizeType i = 0; i < weight_.size(); i++)
					msg<<weight_[i]<<" ";
				progress_.printline(msg,std::cout);

				PsimagLite::OstringStream msg2;
				msg2<<"gsNorm="<<std::norm(commonTargetting_.psi())<<" norms= ";
				for (SizeType i = 0; i < weight_.size(); i++)
					msg2<<normSquared(i)<<" ";
				progress_.printline(msg2,std::cout);
			}

			void printEnergies() const
			{
				for (SizeType i=0;i<commonTargetting_.targetVectors().size();i++)
					printEnergies(commonTargetting_.targetVectors()[i],i);
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
				SizeType threadId = 0;
				typename ModelType::ModelHelperType modelHelper(p,lrs_,threadId);
				typename LanczosSolverType::LanczosMatrixType lanczosHelper(&model_,
				                                                            &modelHelper);


				SizeType total = phi.effectiveSize(i0);
				TargetVectorType phi2(total);
				phi.extract(phi2,i0);
				TargetVectorType x(total);
				lanczosHelper.matrixVectorProduct(x,phi2);
				PsimagLite::OstringStream msg;
				msg<<"Hamiltonian average at time="<<commonTargetting_.currentTime();
				msg<<" for target="<<whatTarget;
				msg<<" sector="<<i0<<" <phi(t)|H|phi(t)>="<<(phi2*x)<<" <phi(t)|phi(t)>="<<(phi2*phi2);
				progress_.printline(msg,std::cout);
			}

			// in situ computation:
			void cocoon(SizeType direction,const BlockType& block) const
			{
				std::cout<<"-------------&*&*&* In-situ measurements start\n";

				if (commonTargetting_.noStageIs(DISABLED))
					std::cout<<"ALL OPERATORS HAVE BEEN APPLIED\n";
				else
					std::cout<<"NOT ALL OPERATORS APPLIED YET\n";

				PsimagLite::String modelName = model_.params().model;

				if (modelName == "HubbardOneBand" ||
				    modelName == "HubbardOneBandExtended" ||
				    modelName == "Immm") {
					commonTargetting_.cocoonLegacy(direction,block);
				}

				SizeType site = block[0];

				commonTargetting_.cocoon(direction,site,commonTargetting_.psi(),"PSI");
				commonTargetting_.cocoon(direction,site,commonTargetting_.targetVectors()[0],"P0");

				std::cout<<"-------------&*&*&* In-situ measurements end\n";
			}

			void checkNorms() const
			{
				PsimagLite::OstringStream msg;
				msg<<"Checking norms: ";
				for (SizeType i=0;i<commonTargetting_.targetVectors().size();i++) {
					RealType norma = std::norm(commonTargetting_.targetVectors()[i]);
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

				TargetVectorType r(n2);

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

			const LeftRightSuperType& lrs_;
			const ModelType& model_;
			const TargettingParamsType& tstStruct_;
			const WaveFunctionTransfType& wft_;
			PsimagLite::ProgressIndicator progress_;
			typename PsimagLite::Vector<RealType>::Type times_,weight_;
			CommonTargettingType commonTargetting_;
			RealType gsWeight_;
	};     //class TimeStepTargetting

	template<template<typename,typename,typename> class LanczosSolverTemplate,
	         typename MatrixVectorType,
	         typename WaveFunctionTransfType,
	         typename IoType_>
	std::ostream& operator<<(std::ostream& os,
	                         const TimeStepTargetting<LanczosSolverTemplate,
	                         MatrixVectorType,
	                         WaveFunctionTransfType,IoType_>& tst)
	{
		tst.print(os);
		return os;
	}
} // namespace Dmrg

#endif
