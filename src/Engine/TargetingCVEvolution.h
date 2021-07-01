/*
Copyright (c) 2009-2014-2021, UT-Battelle, LLC
All rights reserved

[DMRG++, Version 6.]
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

#ifndef TARGETING_CVEvolution_H
#define TARGETING_CVEvolution_H

#include <iostream>
#include "ProgressIndicator.h"
#include "TargetParamsCorrectionVector.h"
#include "ProgramGlobals.h"
#include "ParametersForSolver.h"
#include "TimeVectorsKrylov.h"
#include "TargetingBase.h"
#include "BlockDiagonalMatrix.h"
#include "PredicateAwesome.h"
#include "CorrectionVectorSkeleton.h"

namespace Dmrg {

template<typename LanczosSolverType_, typename VectorWithOffsetType_>
class TargetingCVEvolution : public TargetingBase<LanczosSolverType_,VectorWithOffsetType_> {

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
	typedef TargetParamsCorrectionVector<ModelType> TargetParamsType;
	typedef typename BasisType::BlockType BlockType;
	typedef typename TargetingCommonType::TimeSerializerType TimeSerializerType;
	typedef typename OperatorType::StorageType SparseMatrixType;
	typedef typename ModelType::InputValidatorType InputValidatorType;
	typedef typename BasisType::QnType QnType;
	typedef typename TargetingCommonType::StageEnumType StageEnumType;
	typedef CorrectionVectorSkeleton<LanczosSolverType,
	VectorWithOffsetType,
	BaseType,
	TargetParamsType> CorrectionVectorSkeletonType;

	TargetingCVEvolution(const LeftRightSuperType& lrs,
	                     const ModelType& model,
	                     const WaveFunctionTransfType& wft,
	                     const QnType&,
	                     InputValidatorType& ioIn)
	    : BaseType(lrs, model, wft, 0),
	      tstStruct_(ioIn, "TargetingCVEvolution", model),
	      wft_(wft),
	      progress_("TargetingCVEvolution"),
	      counter_(0),
	      skeleton_(ioIn, tstStruct_, model, lrs, this->common().aoe().energy()),
	      weight_(PsimagLite::IsComplexNumber<ComplexOrRealType>::True ? 3 : 5),
	      gsWeight_(tstStruct_.gsWeight())
	{
		if (!PsimagLite::IsComplexNumber<ComplexOrRealType>::True)
			err("TargetingCVEvolution: only for complex for now (sorry!)\n");

		if (!wft.isEnabled())
			err("TST needs an enabled wft\n");

		if (tstStruct_.sites() == 0)
			err("TST needs at least one TSPSite\n");

		if (gsWeight_ < 0 || gsWeight_ >= 1)
			err("gsWeight_ must be in [0, 1)\n");

		SizeType n = weight_.size();
		RealType factor = (1.0 - gsWeight_)/n;
		for (SizeType i = 0; i < n; ++i) weight_[i] = factor;
	}

	SizeType sites() const { return tstStruct_.sites(); }

	SizeType targets() const { return weight_.size(); }

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

		if (site > 1 && site < numberOfSites - 2)
			return;

		if (direction == ProgramGlobals::DirectionEnum::EXPAND_SYSTEM) {
			if (site == 1) return;
		} else {
			if (site == numberOfSites - 2) return;
		}

		SizeType x = (site == 1) ? 0 : numberOfSites - 1;
		BlockType block(1, x);
		evolveInternal(Eg,direction,block,loopNumber);
	}

	bool end() const
	{
		return (tstStruct_.maxTime() != 0 &&
		        this->common().aoe().time() >= tstStruct_.maxTime());
	}

	void read(typename TargetingCommonType::IoInputType& io, PsimagLite::String prefix)
	{
		this->common().readGSandNGSTs(io, prefix, "TimeStep");
	}

	void write(const VectorSizeType& block,
	           PsimagLite::IoSelector::Out& io,
	           PsimagLite::String prefix) const
	{
		PsimagLite::OstringStream msgg(std::cout.precision());
		PsimagLite::OstringStream::OstringStreamType& msg = msgg();
		msg<<"Saving state...";
		progress_.printline(msgg, std::cout);

		this->common().write(io, block, prefix);
		this->common().writeNGSTs(io, prefix, block, "TimeStep");
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

		if (phiNew.size() == 0) return;

		this->common().aoe().targetVectors(0) = phiNew;
		VectorWithOffsetType bogusTv;

		if (counter_ == 0) {
			skeleton_.calcDynVectors(phiNew,
			                         this->common().aoe().targetVectors(1),
			                         bogusTv);
			skeleton_.calcDynVectors(this->common().aoe().targetVectors(1),
			                         this->common().aoe().targetVectors(2),
			                         bogusTv);
		} else {
			// wft tv1
			skeleton_.calcDynVectors(this->common().aoe().targetVectors(1),
			                         this->common().aoe().targetVectors(2),
			                         bogusTv);
		}

		bool doBorderIfBorder = false;
		this->common().cocoon(block1, direction, doBorderIfBorder);

		this->common().printNormsAndWeights(gsWeight_, weight_);
	}

	TargetParamsType tstStruct_;
	const WaveFunctionTransfType& wft_;
	PsimagLite::ProgressIndicator progress_;
	SizeType counter_;
	CorrectionVectorSkeletonType skeleton_;
	VectorRealType weight_;
	RealType gsWeight_;
};     //class TargetingCVEvolution
} // namespace Dmrg

#endif

