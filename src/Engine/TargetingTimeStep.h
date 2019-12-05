/*
Copyright (c) 2009-2014, UT-Battelle, LLC
All rights reserved

[DMRG++, Version 5.]
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

#ifndef TARGETING_TIMESTEP_H
#define TARGETING_TIMESTEP_H

#include <iostream>
#include "ProgressIndicator.h"
#include "TargetParamsTimeStep.h"
#include "ProgramGlobals.h"
#include "ParametersForSolver.h"
#include "TimeVectorsKrylov.h"
#include "TimeVectorsRungeKutta.h"
#include "TimeVectorsSuzukiTrotter.h"
#include "TargetingBase.h"
#include "BlockDiagonalMatrix.h"
#include "PredicateAwesome.h"

namespace Dmrg {

template<typename LanczosSolverType_, typename VectorWithOffsetType_>
class TargetingTimeStep : public TargetingBase<LanczosSolverType_,VectorWithOffsetType_> {

	enum {BORDER_NEITHER, BORDER_LEFT, BORDER_RIGHT};

public:

	typedef LanczosSolverType_ LanczosSolverType;
	typedef TargetingBase<LanczosSolverType,VectorWithOffsetType_> BaseType;
	typedef typename BaseType::TargetingCommonType TargetingCommonType;
	typedef std::pair<SizeType,SizeType> PairType;
	typedef typename BaseType::OptionsType OptionsType;
	typedef typename BaseType::MatrixVectorType MatrixVectorType;
	typedef typename MatrixVectorType::ModelType ModelType;
	typedef typename ModelType::RealType RealType;
	typedef typename ModelType::OperatorsType OperatorsType;
	typedef typename ModelType::ModelHelperType ModelHelperType;
	typedef typename ModelHelperType::LeftRightSuperType LeftRightSuperType;
	typedef typename LeftRightSuperType::BasisWithOperatorsType BasisWithOperatorsType;
	typedef PsimagLite::Vector<SizeType>::Type VectorSizeType;
	typedef typename BaseType::WaveFunctionTransfType WaveFunctionTransfType;
	typedef typename WaveFunctionTransfType::VectorWithOffsetType VectorWithOffsetType;
	typedef typename VectorWithOffsetType::value_type ComplexOrRealType;
	typedef typename VectorWithOffsetType::VectorType TargetVectorType;
	typedef typename PsimagLite::Vector<RealType>::Type VectorRealType;
	typedef typename BasisWithOperatorsType::OperatorType OperatorType;
	typedef typename BasisWithOperatorsType::BasisType BasisType;
	typedef TargetParamsTimeStep<ModelType> TargetParamsType;
	typedef typename BasisType::BlockType BlockType;
	typedef typename TargetingCommonType::TimeSerializerType TimeSerializerType;
	typedef typename OperatorType::StorageType SparseMatrixType;
	typedef typename ModelType::InputValidatorType InputValidatorType;
	typedef typename BasisType::QnType QnType;
	typedef typename TargetingCommonType::StageEnumType StageEnumType;

	TargetingTimeStep(const LeftRightSuperType& lrs,
	                  const ModelType& model,
	                  const WaveFunctionTransfType& wft,
	                  const QnType&,
	                  InputValidatorType& ioIn,
	                  PsimagLite::String targeting)
	    : BaseType(lrs, model, wft, 0),
	      tstStruct_(ioIn, targeting, model),
	      wft_(wft),
	      progress_(targeting),
	      times_(tstStruct_.timeSteps()),
	      weight_(tstStruct_.timeSteps()),
	      tvEnergy_(times_.size(),0.0),
	      gsWeight_(tstStruct_.gsWeight())
	{
		if (!wft.isEnabled())
			err("TST needs an enabled wft\n");
		if (tstStruct_.sites() == 0)
			err("TST needs at least one TSPSite\n");

		RealType tau =tstStruct_.tau();
		RealType sum = 0;
		SizeType n = times_.size();

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

		this->common().aoe().initTimeVectors(tstStruct_, times_, ioIn);
	}

	SizeType sites() const { return tstStruct_.sites(); }

	SizeType targets() const { return tstStruct_.timeSteps(); }

	RealType weight(SizeType i) const
	{
		assert(!this->common().aoe().allStages(StageEnumType::DISABLED));
		return weight_[i];
	}

	RealType gsWeight() const
	{
		if (this->common().aoe().allStages(StageEnumType::DISABLED))
			return 1.0;
		return gsWeight_;
	}

	bool includeGroundStage() const
	{
		if (!this->common().aoe().noStageIs(StageEnumType::DISABLED))
			return true;
		bool b = (fabs(gsWeight_)>1e-6);
		return b;
	}

	void evolve(const VectorRealType& energies,
	            ProgramGlobals::DirectionEnum direction,
	            const BlockType& block1,
	            const BlockType&,
	            SizeType loopNumber)
	{
		assert(block1.size() > 0);
		SizeType site = block1[0];
		assert(energies.size() > 0);
		RealType Eg = energies[0];
		evolveInternal(Eg,direction,block1,loopNumber);
		SizeType numberOfSites = this->lrs().super().block().size();

		if (site > 1 && site < numberOfSites-2)
			return;
		if (site == 1 && direction == ProgramGlobals::DirectionEnum::EXPAND_SYSTEM)
			return;

		SizeType x = (site == 1) ? 0 : numberOfSites-1;
		BlockType block(1,x);
		evolveInternal(Eg,direction,block,loopNumber);
	}

	bool end() const
	{
		return (tstStruct_.maxTime() != 0 &&
		        this->common().aoe().time() >= tstStruct_.maxTime());
	}

	void read(typename TargetingCommonType::IoInputType& io, PsimagLite::String prefix)
	{
		this->common().readGSandNGSTs(io, prefix);
	}

	void write(const VectorSizeType& block,
	           PsimagLite::IoSelector::Out& io,
	           PsimagLite::String prefix) const
	{
		PsimagLite::OstringStream msg;
		msg<<"Saving state...";
		progress_.printline(msg,std::cout);

		this->common().write(io, block, prefix);
		this->common().writeNGSTs(io, block, prefix);
	}

private:

	void evolveInternal(RealType Eg,
	                    ProgramGlobals::DirectionEnum direction,
	                    const BlockType& block1,
	                    SizeType loopNumber)
	{
		if (direction == ProgramGlobals::DirectionEnum::INFINITE) return;
		VectorWithOffsetType phiNew;
		assert(block1.size() > 0);
		SizeType site = block1[0];
		this->common().aoe().getPhi(&phiNew,
		                            Eg,
		                            direction,
		                            site,
		                            loopNumber,
		                            tstStruct_);

		PairType startEnd(0,times_.size());
		bool allOperatorsApplied = (this->common().aoe().noStageIs(StageEnumType::DISABLED) &&
		                            this->common().aoe().noStageIs(StageEnumType::OPERATOR));

		VectorSizeType indices(startEnd.second - startEnd.first);
		for (SizeType i = 0; i < indices.size(); ++i) indices[i] = i + startEnd.first;

		const bool isLastCall = true;
		this->common().aoe().calcTimeVectors(indices,
		                                     Eg,
		                                     phiNew,
		                                     direction,
		                                     allOperatorsApplied,
		                                     false, // don't wft or advance indices[0]
		                                     block1,
		                                     isLastCall);

		bool doBorderIfBorder = false;
		this->common().cocoon(block1, direction, doBorderIfBorder);

		PsimagLite::String predicate = this->model().params().printHamiltonianAverage;
		const SizeType center = this->model().geometry().numberOfSites()/2;
		PsimagLite::PredicateAwesome<>::replaceAll(predicate, "c", ttos(center));
		PsimagLite::PredicateAwesome<> pAwesome(predicate);
		assert(block1.size() > 0);
		if (pAwesome.isTrue("s", block1[0]))
			printEnergies(); // in-situ

		const OptionsType& options = this->model().params().options;
		bool normalizeTimeVectors = (options.isSet("normalizeTimeVectors") ||
		                             options.isSet("TargetingAncilla"));

		if (options.isSet("neverNormalizeVectors"))
			normalizeTimeVectors = false;

		if (normalizeTimeVectors)
			this->common().normalizeTimeVectors();

		this->common().printNormsAndWeights(gsWeight_, weight_);
	}

	void printEnergies() const
	{
		for (SizeType i=0;i<this->common().aoe().targetVectors().size();i++)
			printEnergies(this->common().aoe().targetVectors()[i],i);
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
		SizeType p = this->lrs().super().findPartitionNumber(phi.offset(i0));
		typename ModelHelperType::Aux aux(p, BaseType::lrs());
		typename ModelType::HamiltonianConnectionType hc(BaseType::lrs(),
		                                                 BaseType::model().geometry(),
		                                                 ModelType::modelLinks(),
		                                                 this->common().aoe().time(),
		                                                 0);
		typename LanczosSolverType::MatrixType lanczosHelper(BaseType::model(),
		                                                     hc,
		                                                     aux);

		SizeType total = phi.effectiveSize(i0);
		TargetVectorType phi2(total);
		phi.extract(phi2,i0);
		TargetVectorType x(total);
		lanczosHelper.matrixVectorProduct(x,phi2);
		PsimagLite::OstringStream msg;
		msg<<"Hamiltonian average at time="<<this->common().aoe().time();
		msg<<" for target="<<whatTarget;
		ComplexOrRealType numerator = phi2*x;
		ComplexOrRealType den = phi2*phi2;
		ComplexOrRealType division = (PsimagLite::norm(den)<1e-10) ? 0 : numerator/den;
		msg<<" sector="<<i0<<" <phi(t)|H|phi(t)>="<<numerator;
		msg<<" <phi(t)|phi(t)>="<<den<<" "<<division;
		progress_.printline(msg,std::cout);
		tvEnergy_[whatTarget] = PsimagLite::real(division);
	}

	TargetParamsType tstStruct_;
	const WaveFunctionTransfType& wft_;
	PsimagLite::ProgressIndicator progress_;
	VectorRealType times_;
	VectorRealType weight_;
	mutable VectorRealType tvEnergy_;
	RealType gsWeight_;
};     //class TargetingTimeStep
} // namespace Dmrg

#endif

