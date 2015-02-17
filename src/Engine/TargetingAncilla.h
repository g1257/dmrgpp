/*
Copyright (c) 2009-2015, UT-Battelle, LLC
All rights reserved

[DMRG++, Version 3.0]
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

#ifndef DMRG_TARGETING_ANCILLA_H
#define DMRG_TARGETING_ANCILLA_H

#include <iostream>
#include "ProgressIndicator.h"
#include "BLAS.h"
#include "TimeSerializer.h"
#include "TargetParamsTimeStep.h"
#include "ProgramGlobals.h"
#include "ParametersForSolver.h"
#include "ParallelWft.h"
#include "TimeVectorsKrylov.h"
#include "TimeVectorsRungeKutta.h"
#include "TimeVectorsSuzukiTrotter.h"
#include "TargetingBase.h"
#include "BlockMatrix.h"

namespace Dmrg {

template<template<typename,typename,typename> class LanczosSolverTemplate,
         typename MatrixVectorType_,
         typename WaveFunctionTransfType_>
class TargetingAncilla : public TargetingBase<LanczosSolverTemplate,
                                               MatrixVectorType_,
                                               WaveFunctionTransfType_> {

	enum {BORDER_NEITHER, BORDER_LEFT, BORDER_RIGHT};

public:

	typedef TargetingBase<LanczosSolverTemplate,
	                      MatrixVectorType_,
	                      WaveFunctionTransfType_> BaseType;
	typedef std::pair<SizeType,SizeType> PairType;
	typedef MatrixVectorType_ MatrixVectorType;
	typedef typename MatrixVectorType::ModelType ModelType;
	typedef typename ModelType::RealType RealType;
	typedef std::complex<RealType> ComplexType;
	typedef typename ModelType::OperatorsType OperatorsType;
	typedef typename ModelType::ModelHelperType ModelHelperType;
	typedef typename ModelHelperType::LeftRightSuperType LeftRightSuperType;
	typedef typename LeftRightSuperType::BasisWithOperatorsType BasisWithOperatorsType;
	typedef PsimagLite::Vector<SizeType>::Type VectorSizeType;
	typedef WaveFunctionTransfType_ WaveFunctionTransfType;
	typedef typename WaveFunctionTransfType::VectorWithOffsetType VectorWithOffsetType;
	typedef typename VectorWithOffsetType::VectorType TargetVectorType;
	typedef PsimagLite::ParametersForSolver<RealType> ParametersForSolverType;
	typedef LanczosSolverTemplate<ParametersForSolverType,
	                              MatrixVectorType,
	                              TargetVectorType> LanczosSolverType;
	typedef typename PsimagLite::Vector<RealType>::Type VectorType;
	typedef PsimagLite::Matrix<ComplexType> ComplexMatrixType;
	typedef typename BasisWithOperatorsType::OperatorType OperatorType;
	typedef typename BasisWithOperatorsType::BasisType BasisType;
	typedef TargetParamsTimeStep<ModelType> TargettingParamsType;
	typedef typename BasisType::BlockType BlockType;
	typedef BlockMatrix<ComplexMatrixType> ComplexBlockMatrixType;
	typedef TimeSerializer<VectorWithOffsetType> TimeSerializerType;
	typedef typename OperatorType::SparseMatrixType SparseMatrixType;
	typedef typename BasisWithOperatorsType::BasisDataType BasisDataType;
	typedef typename ModelType::InputValidatorType InputValidatorType;

	enum {DISABLED,OPERATOR,WFT_NOADVANCE,WFT_ADVANCE};

	static SizeType const PRODUCT = TargettingParamsType::PRODUCT;
	static SizeType const SUM = TargettingParamsType::SUM;

	TargetingAncilla(const LeftRightSuperType& lrs,
	                  const ModelType& model,
	                  const TargettingParamsType& tstStruct,
	                  const WaveFunctionTransfType& wft,
	                  const SizeType&,
	                  InputValidatorType& ioIn)
	    : BaseType(lrs,model,tstStruct,wft,tstStruct.timeSteps(),0),
	      tstStruct_(tstStruct),
	      wft_(wft),
	      progress_("TargetingAncilla"),
	      times_(tstStruct_.timeSteps()),
	      weight_(tstStruct_.timeSteps())
	{
		if (!wft.isEnabled()) {
			PsimagLite::String msg("TargetingAncilla needs an enabled wft\n");
			throw PsimagLite::RuntimeError(msg);
		}

		if (tstStruct_.sites() == 0) {
			PsimagLite::String msg("TargetingAncilla needs at least one TSPSite\n");
			throw PsimagLite::RuntimeError(msg);
		}

		RealType tau =tstStruct_.tau();
		RealType sum = 0;
		SizeType n = times_.size();
		gsWeight_ = (tstStruct_.concatenation() == SUM) ? 0.1 : 0.0;
		gsWeight_ = this->common().setGsWeight(gsWeight_);

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

		this->common().initTimeVectors(times_,ioIn);
	}

	RealType weight(SizeType i) const
	{
		assert(!this->common().allStages(DISABLED));
		return weight_[i];
	}

	RealType gsWeight() const
	{
		if (this->common().allStages(DISABLED)) return 1.0;
		return gsWeight_;
	}

	bool includeGroundStage() const
	{
		if (!this->common().noStageIs(DISABLED)) return true;
		bool b = (fabs(gsWeight_)>1e-6);
		return b;
	}

	void evolve(RealType Eg,
	            SizeType direction,
	            const BlockType& block1,
	            const BlockType&,
	            SizeType loopNumber)
	{
		VectorWithOffsetType phiNew;
		this->common().getPhi(phiNew,Eg,direction,block1[0],loopNumber);

		PairType startEnd(0,times_.size());
		bool allOperatorsApplied = (this->common().noStageIs(DISABLED) &&
		                            this->common().noStageIs(OPERATOR));

		this->common().calcTimeVectors(startEnd,
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
		this->common().template load<TimeSerializerType>(f);
	}

	void print(std::ostream& os) const
	{
		os<<"TSTWeightsTimeVectors=";
		for (SizeType i=0;i<weight_.size();i++)
			os<<weight_[i]<<" ";
		os<<"\n";
		os<<"TSTWeightGroundState="<<gsWeight_<<"\n";
	}

	template<typename IoOutputType>
	void save(const VectorSizeType& block,IoOutputType& io) const
	{
		PsimagLite::OstringStream msg;
		msg<<"Saving state...";
		progress_.printline(msg,std::cout);

		SizeType marker = 0;
		if (this->common().noStageIs(DISABLED)) marker = 1;

		TimeSerializerType ts(this->common().currentTime(),
		                      block[0],
		        this->common().targetVectors(),
		        marker);
		ts.save(io);
		this->common().psi().save(io,"PSI");
	}

	void updateOnSiteForTimeDep(BasisWithOperatorsType& basisWithOps) const
	{

		BlockType X = basisWithOps.block();
		if (X.size()!=1) return;
		assert(X[0]==0 || X[0]==this->leftRightSuper().super().block().size()-1);
		typename PsimagLite::Vector<OperatorType>::Type creationMatrix;
		SparseMatrixType hmatrix;
		BasisDataType q;
		this->model().setNaturalBasis(creationMatrix,
		                              hmatrix,
		                              q,
		                              X,
		                              this->common().currentTime());
		basisWithOps.setVarious(X,hmatrix,q,creationMatrix);
	}

	bool end() const
	{
		return (tstStruct_.maxTime() != 0 &&
		        this->common().currentTime() >= tstStruct_.maxTime());
	}

private:

	void printNormsAndWeights() const
	{
		if (this->common().allStages(DISABLED)) return;

		PsimagLite::OstringStream msg;
		msg<<"gsWeight="<<gsWeight_<<" weights= ";
		for (SizeType i = 0; i < weight_.size(); i++)
			msg<<weight_[i]<<" ";
		progress_.printline(msg,std::cout);

		PsimagLite::OstringStream msg2;
		msg2<<"gsNorm="<<std::norm(this->common().psi())<<" norms= ";
		for (SizeType i = 0; i < weight_.size(); i++)
			msg2<<this->common().normSquared(i)<<" ";
		progress_.printline(msg2,std::cout);
	}

	void printEnergies() const
	{
		for (SizeType i=0;i<this->common().targetVectors().size();i++)
			printEnergies(this->common().targetVectors()[i],i);
	}

	void printEnergies(const VectorWithOffsetType& phi,SizeType whatTarget) const
	{
		for (SizeType ii=0;ii<phi.sectors();ii++) {
			SizeType i = phi.sector(ii);
			printEnergies(phi,whatTarget,i);
		}
	}

	void printEnergies(const VectorWithOffsetType& phi,
	                   SizeType whatTarget,
	                   SizeType i0) const
	{
		SizeType p = this->leftRightSuper().super().findPartitionNumber(phi.offset(i0));
		SizeType threadId = 0;
		typename ModelType::ModelHelperType modelHelper(p,
		                                                this->leftRightSuper(),
		                                                threadId);
		typename LanczosSolverType::LanczosMatrixType lanczosHelper(&this->model(),
		                                                            &modelHelper);

		SizeType total = phi.effectiveSize(i0);
		TargetVectorType phi2(total);
		phi.extract(phi2,i0);
		TargetVectorType x(total);
		lanczosHelper.matrixVectorProduct(x,phi2);
		PsimagLite::OstringStream msg;
		msg<<"Hamiltonian average at time="<<this->common().currentTime();
		msg<<" for target="<<whatTarget;
		msg<<" sector="<<i0<<" <phi(t)|H|phi(t)>="<<(phi2*x);
		msg<<" <phi(t)|phi(t)>="<<(phi2*phi2);
		progress_.printline(msg,std::cout);
	}

	// in situ computation:
	void cocoon(SizeType direction,const BlockType& block) const
	{
		std::cout<<"-------------&*&*&* In-situ measurements start\n";

		if (this->common().noStageIs(DISABLED))
			std::cout<<"ALL OPERATORS HAVE BEEN APPLIED\n";
		else
			std::cout<<"NOT ALL OPERATORS APPLIED YET\n";

		PsimagLite::String modelName = this->model().params().model;

		if (modelName == "HubbardOneBand" ||
		    modelName == "HubbardOneBandExtended" ||
		    modelName == "Immm") {
			this->common().cocoonLegacy(direction,block);
		}

		this->common().cocoon(block,direction);

		std::cout<<"-------------&*&*&* In-situ measurements end\n";
	}

	const TargettingParamsType& tstStruct_;
	const WaveFunctionTransfType& wft_;
	PsimagLite::ProgressIndicator progress_;
	typename PsimagLite::Vector<RealType>::Type times_,weight_;
	RealType gsWeight_;
};     //class TargetingAncilla

template<template<typename,typename,typename> class LanczosSolverTemplate,
         typename MatrixVectorType,
         typename WaveFunctionTransfType>
std::ostream& operator<<(std::ostream& os,
                         const TargetingAncilla<LanczosSolverTemplate,
                         MatrixVectorType,
                         WaveFunctionTransfType>& tst)
{
	tst.print(os);
	return os;
}
} // namespace Dmrg

#endif

