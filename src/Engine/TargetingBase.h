/*
Copyright (c) 2014, UT-Battelle, LLC
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

/*! \file TargetingBase.h
 *
 *
 *
 */

#ifndef TARGETING_BASE_H
#define TARGETING_BASE_H
#include <iostream>
#include "TargetParamsBase.h"
#include "TargetHelper.h"
#include "TargetingCommon.h"
#include "Wft/WaveFunctionTransfFactory.h"
#include "IoSelector.h"

namespace Dmrg {

template<typename LanczosSolverType_, typename VectorWithOffsetType_>
class TargetingBase {

public:

	typedef LanczosSolverType_ LanczosSolverType;
	typedef VectorWithOffsetType_ VectorWithOffsetType;
	typedef typename LanczosSolverType::LanczosMatrixType MatrixVectorType;
	typedef typename MatrixVectorType::ModelType ModelType;
	typedef typename ModelType::RealType RealType;
	typedef PsimagLite::ParametersForSolver<RealType> ParametersForSolverType;
	typedef typename ModelType::ModelHelperType ModelHelperType;
	typedef typename ModelHelperType::LeftRightSuperType LeftRightSuperType;
	typedef typename LeftRightSuperType::BasisWithOperatorsType BasisWithOperatorsType;
	typedef typename BasisWithOperatorsType::SparseMatrixType SparseMatrixType;
	typedef typename BasisWithOperatorsType::OperatorType OperatorType;
	typedef typename BasisWithOperatorsType::BasisType BasisType;
	typedef typename SparseMatrixType::value_type ComplexOrRealType;
	typedef typename BasisType::BlockType BlockType;
	typedef WaveFunctionTransfFactory<LeftRightSuperType,
	                                  VectorWithOffsetType> WaveFunctionTransfType;
	typedef typename VectorWithOffsetType::VectorType VectorType;
	typedef VectorType TargetVectorType;
	typedef TargetParamsBase<ModelType> TargetParamsType;
	typedef TargetHelper<ModelType,
	                     TargetParamsType,
	                     WaveFunctionTransfType> TargetHelperType;
	typedef TargetingCommon<TargetHelperType,
	                        VectorWithOffsetType,
	                        LanczosSolverType> TargetingCommonType;
	typedef typename TargetingCommonType::ApplyOperatorExpressionType
	ApplyOperatorExpressionType;
	typedef typename BasisWithOperatorsType::SymmetryElectronsSzType
	SymmetryElectronsSzType;
	typedef typename PsimagLite::Vector<OperatorType>::Type VectorOperatorType;

	enum {DISABLED=ApplyOperatorExpressionType::DISABLED,
		  OPERATOR=ApplyOperatorExpressionType::OPERATOR,
		  WFT_NOADVANCE=ApplyOperatorExpressionType::WFT_NOADVANCE};

	TargetingBase(const LeftRightSuperType& lrs,
	              const ModelType& model,
	              const WaveFunctionTransfType& wft,
	              SizeType indexNoAdvance)
	    : lrs_(lrs),
	      model_(model),
	      commonTargeting_(lrs,model,wft,indexNoAdvance)
	{}

	virtual ~TargetingBase() {}

	virtual RealType gsWeight() const = 0;

	virtual RealType weight(SizeType i) const = 0;

	virtual void evolve(RealType Eg,
	                    ProgramGlobals::DirectionEnum direction,
	                    const BlockType& block1,
	                    const BlockType& block2,
	                    SizeType loopNumber) = 0;

	virtual bool includeGroundStage() const {return true; }

	virtual void updateOnSiteForCorners(BasisWithOperatorsType& basisWithOps) const
	{
		if (BasisWithOperatorsType::useSu2Symmetry()) return;

		BlockType X = basisWithOps.block();

		if (X.size()!=1) return;

		if (X[0] != 0 && X[0] != lrs_.super().block().size()-1)
			return;

		basisWithOps.setVarious(X, model_, commonTargeting_.currentTime());
	}

	virtual bool end() const
	{
		return false;
	}

	virtual SizeType size() const
	{
		if (commonTargeting_.allStages(DISABLED)) return 0;
		return commonTargeting_.targetVectors().size();
	}

	virtual RealType normSquared(SizeType i) const
	{
		return commonTargeting_.normSquared(i);
	}

	virtual void load(const PsimagLite::String&) = 0;

	virtual void save(const typename PsimagLite::Vector<SizeType>::Type&,
	                  PsimagLite::IoSelector::Out&) const = 0;

	// non-virtual below

	const ModelType& model() const { return model_; }

	template<typename SomeBasisType>
	void setGs(const typename PsimagLite::Vector<VectorType>::Type& v,
	           const SomeBasisType& someBasis)
	{
		commonTargeting_.setGs(v,someBasis);
	}

	const VectorWithOffsetType& gs() const
	{
		return commonTargeting_.psi();
	}

	const VectorWithOffsetType& operator()(SizeType i) const
	{
		return commonTargeting_.targetVectors()[i];
	}

	void initialGuess(VectorWithOffsetType& initialVector,
	                  const typename PsimagLite::Vector<SizeType>::Type& block,
	                  bool noguess) const
	{
		commonTargeting_.initialGuess(initialVector, block, noguess);
	}

	const RealType& time() const {return commonTargeting_.currentTime(); }

	const ComplexOrRealType& inSitu(SizeType i) const
	{
		return commonTargeting_.inSitu(i);
	}

	const LeftRightSuperType& lrs() const { return lrs_; }

	template<typename IoOutType>
	void print(IoOutType&) const
	{
		std::cerr<<__FILE__<<" print() isn't needed WARNING REMOVE FIXME\n";
	}

protected:

	TargetingCommonType& common()
	{
		return commonTargeting_;
	}

	const TargetingCommonType& common() const
	{
		return commonTargeting_;
	}

private:

	const LeftRightSuperType& lrs_;
	const ModelType& model_;
	TargetingCommonType commonTargeting_;
};     //class TargetingBase

} // namespace Dmrg
/*@}*/
#endif

