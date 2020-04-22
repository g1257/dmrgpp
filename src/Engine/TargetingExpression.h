/*
Copyright (c) 2009-2013-2019, UT-Battelle, LLC
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

/** \ingroup DMRG */
/*@{*/

/*! \file TargetingExpression.h
 * TBW FIXME TODO
 */

#ifndef TARGETING_EXPRESSION_H
#define TARGETING_EXPRESSION_H
#include <iostream>
#include "TargetingBase.h"
#include <stdexcept>
#include "Pvector.h"
#include "SpecForTargetingExpression.h"
#include "CanonicalExpression.h"

namespace Dmrg {

template<typename LanczosSolverType_, typename VectorWithOffsetType_>
class TargetingExpression : public TargetingBase<LanczosSolverType_,VectorWithOffsetType_> {

	typedef TargetingBase<LanczosSolverType_,VectorWithOffsetType_> BaseType;
	typedef typename BaseType::TargetingCommonType TargetingCommonType;
	typedef typename BaseType::WaveFunctionTransfType WaveFunctionTransfType;
	typedef typename BaseType::ModelType ModelType;
	typedef typename ModelType::RealType RealType;
	typedef typename ModelType::InputValidatorType InputValidatorType;
	typedef typename ModelType::ModelHelperType ModelHelperType;
	typedef typename ModelHelperType::LeftRightSuperType LeftRightSuperType;
	typedef typename LeftRightSuperType::BasisWithOperatorsType BasisWithOperatorsType;
	typedef typename BasisWithOperatorsType::BasisType BasisType;
	typedef typename BasisType::BlockType BlockType;
	typedef typename BasisType::QnType QnType;
	typedef typename TargetingCommonType::StageEnumType StageEnumType;
	typedef Pvector<VectorWithOffsetType_> PvectorType;
	typedef typename PsimagLite::Vector<PvectorType*>::Type VectorPvectorType;
	typedef SpecForTargetingExpression<BaseType> SpecForTargetingExpressionType;
	typedef typename SpecForTargetingExpressionType::AlgebraType AlgebraType;
	typedef PsimagLite::CanonicalExpression<SpecForTargetingExpressionType>
	CanonicalExpressionType;
	typedef AuxForTargetingExpression<BaseType> AuxForTargetingExpressionType;
	typedef typename TargetingCommonType::VectorRealType VectorRealType;

public:

	TargetingExpression(const LeftRightSuperType& lrs,
	                    const ModelType& model,
	                    const WaveFunctionTransfType& wft,
	                    const QnType&,
	                    InputValidatorType& io)
	    : BaseType(lrs,model,wft,0),
	      progress_("TargetingExpression"),
	      gsWeight_(0.3)
	{
		io.readline(gsWeight_, "GsWeight=");
		pvectorsFromInput(io);
	}

	~TargetingExpression()
	{
		for (SizeType i = 0; i < pVectors_.size(); ++i) {
			delete pVectors_[i];
			pVectors_[i] = 0;
		}
	}

	SizeType sites() const { return 0; }

	SizeType targets() const { return pVectors_.size(); }

	RealType normSquared(SizeType i) const
	{
		return PsimagLite::real(this->common().aoe().targetVectors()[i]*
		                        this->common().aoe().targetVectors()[i]);
	}

	RealType weight(SizeType i) const
	{
		assert(this->common().aoe().noStageIs(StageEnumType::DISABLED));
		assert(i < pVectors_.size());
		return pVectors_[i]->weight();
	}

	RealType gsWeight() const
	{
		return (this->common().aoe().noStageIs(StageEnumType::DISABLED)) ? gsWeight_ : 1.0;
	}

	void evolve(const VectorRealType&,
	            ProgramGlobals::DirectionEnum direction,
	            const BlockType& block1,
	            const BlockType&,
	            SizeType)
	{
		if (direction == ProgramGlobals::DirectionEnum::INFINITE) return;

		this->common().setAllStagesTo(StageEnumType::WFT_NOADVANCE);
		computePvectors(direction);

		SizeType n = pVectors_.size();
		VectorRealType weight(n);
		for (SizeType i = 0; i < n; ++i)
			weight[i] = pVectors_[i]->weight();

		this->common().printNormsAndWeights(gsWeight_, weight);

		bool doBorderIfBorder = true;
		this->common().cocoon(block1, direction, doBorderIfBorder); // in-situ
	}

	void read(typename TargetingCommonType::IoInputType& io,
	          PsimagLite::String prefix)
	{
		this->common().read(io, prefix);
	}

	void write(const typename PsimagLite::Vector<SizeType>::Type& block,
	           PsimagLite::IoSelector::Out& io,
	           PsimagLite::String prefix) const
	{
		this->common().write(io, block, prefix);
	}

private:

	void pvectorsFromInput(InputValidatorType& io)
	{
		SizeType total = 0;
		io.readline(total, "Pvectors=");
		pVectors_.resize(total);
		PsimagLite::String tmp;
		RealType sum = 0.0;
		for (SizeType i = 0; i < total; ++i) {
			io.readline(tmp, "P" + ttos(i) + "=");
			pVectors_[i] = new PvectorType(tmp);
			sum += pVectors_[i]->weight();
		}

		if (sum == 0.0) return;

		RealType factor = (1.0 - gsWeight_)/sum;
		for (SizeType i = 0; i < total; ++i)
			pVectors_[i]->multiplyWeight(factor);
	}

	void computePvectors(ProgramGlobals::DirectionEnum dir)
	{
		CanonicalExpressionType canonicalExpression(opSpec_);
		SizeType total = pVectors_.size();

		AuxForTargetingExpressionType aux(this->common().aoe(),
		                                  this->model(),
		                                  this->lrs(),
		                                  this->common().aoe().psiConst(),
		                                  this->common().aoe().targetVectors(),
		                                  dir);
		const AlgebraType opEmpty(aux);
		for (SizeType i = 0; i < total; ++i) {
			AlgebraType tmp(aux);
			canonicalExpression(tmp, pVectors_[i]->toString(), opEmpty, aux);
			VectorWithOffsetType_& dst = this->common().aoe().targetVectors(i);
			tmp.finalize(&dst);
			pVectors_[i]->setString(tmp.string());
		}
	}

	PsimagLite::ProgressIndicator progress_;
	RealType gsWeight_;
	VectorPvectorType pVectors_;
	SpecForTargetingExpressionType opSpec_;
};     //class TargetingExpression
} // namespace Dmrg
/*@}*/
#endif

