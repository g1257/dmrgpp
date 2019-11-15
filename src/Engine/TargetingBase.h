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
#include "Io/IoSelector.h"
#include "Intent.h"

namespace Dmrg {

template<typename LanczosSolverType_, typename VectorWithOffsetType_>
class TargetingBase {

public:

	typedef LanczosSolverType_ LanczosSolverType;
	typedef VectorWithOffsetType_ VectorWithOffsetType;
	typedef typename LanczosSolverType::MatrixType MatrixVectorType;
	typedef typename MatrixVectorType::ModelType ModelType;
	typedef typename ModelType::RealType RealType;
	typedef typename ModelType::ParametersType ParametersType;
	typedef typename ParametersType::OptionsType OptionsType;
	typedef PsimagLite::ParametersForSolver<RealType> ParametersForSolverType;
	typedef typename ModelType::ModelHelperType ModelHelperType;
	typedef typename ModelHelperType::LeftRightSuperType LeftRightSuperType;
	typedef typename LeftRightSuperType::BasisWithOperatorsType BasisWithOperatorsType;
	typedef typename BasisWithOperatorsType::SparseMatrixType SparseMatrixType;
	typedef typename BasisWithOperatorsType::OperatorType OperatorType;
	typedef typename BasisWithOperatorsType::BasisType BasisType;
	typedef typename SparseMatrixType::value_type ComplexOrRealType;
	typedef typename BasisType::BlockType BlockType;
	typedef typename BasisType::QnType QnType;
	typedef WaveFunctionTransfFactory<LeftRightSuperType, VectorWithOffsetType, OptionsType>
	WaveFunctionTransfType;
	typedef typename VectorWithOffsetType::VectorType VectorType;
	typedef VectorType TargetVectorType;
	typedef TargetParamsBase<ModelType> TargetParamsType;
	typedef TargetHelper<ModelType, WaveFunctionTransfType> TargetHelperType;
	typedef TargetingCommon<TargetHelperType,
	VectorWithOffsetType,
	LanczosSolverType> TargetingCommonType;
	typedef typename TargetingCommonType::ApplyOperatorExpressionType ApplyOperatorExpressionType;
	typedef typename PsimagLite::Vector<OperatorType>::Type VectorOperatorType;
	typedef typename ApplyOperatorExpressionType::StageEnumType StageEnumType;
	typedef typename ApplyOperatorExpressionType::DmrgSerializerType DmrgSerializerType;
	typedef typename PsimagLite::Vector<SizeType>::Type VectorSizeType;
	typedef typename PsimagLite::Vector<RealType>::Type VectorRealType;
	typedef typename PsimagLite::Vector<TargetVectorType>::Type VectorVectorType;
	typedef typename PsimagLite::Vector<VectorVectorType>::Type VectorVectorVectorType;
	typedef typename TargetingCommonType::VectorVectorVectorWithOffsetType
	VectorVectorVectorWithOffsetType;

	TargetingBase(const LeftRightSuperType& lrs,
	              const ModelType& model,
	              const WaveFunctionTransfType& wft,
	              SizeType indexNoAdvance)
	    : lrs_(lrs),
	      model_(model),
	      commonTargeting_(lrs,model,wft,indexNoAdvance)
	{
		Intent<ModelType> intent(model_);
		intent.check();

		const SizeType nexcited = model_.params().numberOfExcited;

		if (nexcited == 1) return; // EARLY EXIT

		PsimagLite::String msg = "nexcited = " + ttos(nexcited) + " > 1 is experimental\n";
		std::cerr<<msg;
		std::cout<<msg;
	}

	TargetingBase(const TargetingBase&) = delete;

	TargetingBase& operator=(const TargetingBase&) = delete;

	virtual ~TargetingBase() {}

	virtual void postCtor()
	{
		commonTargeting_.postCtor(sites(), targets());
	}

	virtual SizeType sites() const = 0;

	virtual SizeType targets() const = 0;

	virtual RealType gsWeight() const = 0;

	virtual RealType weight(SizeType i) const = 0;

	virtual void evolve(const VectorRealType& energies,
	                    ProgramGlobals::DirectionEnum direction,
	                    const BlockType& block1,
	                    const BlockType& block2,
	                    SizeType loopNumber) = 0;

	virtual void read(typename TargetingCommonType::IoInputType&,
	                  PsimagLite::String) = 0;

	virtual void write(const VectorSizeType&,
	                   PsimagLite::IoSelector::Out&,
	                   PsimagLite::String) const = 0;

	// virtuals with default implementation

	virtual bool includeGroundStage() const { return true; }

	virtual void set(VectorVectorVectorType& inV,
	                 const VectorSizeType& sectors,
	                 const BasisType& someBasis)
	{
		const SizeType nsectors = sectors.size();
		const SizeType nexcited = model_.params().numberOfExcited;

		assert(nsectors > 0);

		assert(nexcited > 0);

		if (nsectors != inV.size())
			err("FATAL: inV.size == " + ttos(inV.size()) + " but params.excited.size= "
			    + ttos(nsectors) + "\n");

		for (SizeType sectorIndex = 0; sectorIndex < nsectors; ++sectorIndex) {
			if (inV[sectorIndex].size() != nexcited)
				err("Expected inV[" + ttos(sectorIndex) + "].size == " +
				    ttos(nexcited) + " but found " + ttos(inV[sectorIndex].size()) +
				    " instead\n");

			for (SizeType excitedIndex = 0; excitedIndex < nexcited; ++excitedIndex) {

				commonTargeting_.aoe().setPsi(sectorIndex,
				                              excitedIndex,
				                              inV[sectorIndex][excitedIndex],
				                              someBasis,
				                              sectors);
			}
		}
	}

	virtual void updateOnSiteForCorners(BasisWithOperatorsType& basisWithOps) const
	{
		if (BasisWithOperatorsType::useSu2Symmetry()) return;

		BlockType X = basisWithOps.block();

		if (X.size()!=1) return;

		if (X[0] != 0 && X[0] != lrs_.super().block().size()-1)
			return;

		basisWithOps.setVarious(X, model_, commonTargeting_.aoe().time());
	}

	virtual bool end() const
	{
		return false;
	}

	virtual SizeType size() const
	{
		if (commonTargeting_.aoe().allStages(StageEnumType::DISABLED)) return 0;
		return commonTargeting_.aoe().targetVectors().size();
	}

	virtual RealType normSquared(SizeType i) const
	{
		return commonTargeting_.normSquared(i);
	}

	virtual void initPsi(SizeType nsectors, SizeType nexcited)
	{
		commonTargeting_.aoe().initPsi(nsectors, nexcited);
	}

	// legacy thing for vectorwithoffsets
	virtual void initialGuess(VectorVectorType& initialVector,
	                          const VectorSizeType& block,
	                          bool noguess,
	                          const VectorSizeType& compactedWeights,
	                          const VectorSizeType& sectors,
	                          const BasisType& basis) const
	{
		if (VectorWithOffsetType::name() != "vectorwithoffsets")
			err("FATAL: Wrong execution path\n");

		const VectorWithOffsetType& psi00 = commonTargeting_.aoe().
		        ensureOnlyOnePsi("initialGuess");
		VectorWithOffsetType vwo(compactedWeights, sectors, basis);
		commonTargeting_.initialGuess(vwo, psi00, block, noguess);
		const SizeType n = vwo.sectors();
		initialVector.resize(n);
		for (SizeType i = 0; i < n; ++i)
			vwo.extract(initialVector[i], vwo.sector(i));
	}

	virtual void initialGuess(VectorType& initialVector,
	                          const VectorSizeType& block,
	                          bool noguess,
	                          const VectorSizeType& compactedWeights,
	                          const VectorSizeType& sectors,
	                          SizeType sectorIndex,
	                          SizeType excited,
	                          const BasisType& basis) const
	{
		if (VectorWithOffsetType::name() == "vectorwithoffsets")
			err("FATAL: Wrong execution path\n");

		const VectorVectorVectorWithOffsetType& psi = commonTargeting_.aoe().psiConst();
		const SizeType nsectors = psi.size();

		if (nsectors != compactedWeights.size())
			err("initialGuess compactedWeights\n");

		const SizeType numberOfExcited = psi[sectorIndex].size();

		if (excited > numberOfExcited)
			err("initialGuess, excited=" + ttos(excited) + " > " +
			    ttos(numberOfExcited) + "\n");

		SizeType start = 0;
		SizeType end = numberOfExcited;
		if (excited < numberOfExcited) {
			start = excited;
			end = excited + 1;
		}

		for (SizeType e = start; e < end; ++e) {

			VectorWithOffsetType vwo(compactedWeights[sectorIndex],
			                         sectors[sectorIndex],
			                         basis);
			commonTargeting_.initialGuess(vwo,
			                              *(psi[sectorIndex][e]),
			                              block,
			                              noguess);

			VectorType tmpVector;
			vwo.extract(tmpVector, vwo.sector(0));
			if (e == start)
				initialVector = tmpVector;
			else
				initialVector += tmpVector;
		}
	}

	// non-virtual below

	const ModelType& model() const { return model_; }

	const VectorVectorVectorWithOffsetType& psiConst() const
	{
		return commonTargeting_.aoe().psiConst();
	}

	const VectorWithOffsetType& operator()(SizeType i) const
	{
		return commonTargeting_.aoe().targetVectors()[i];
	}

	RealType time() const {return commonTargeting_.aoe().time(); }

	const ComplexOrRealType& inSitu(SizeType i) const
	{
		return commonTargeting_.inSitu(i);
	}

	const LeftRightSuperType& lrs() const { return lrs_; }

	static PsimagLite::String buildPrefix(PsimagLite::IoSelector::Out& io,
	                                      SizeType counter)
	{
		PsimagLite::String prefix("TargetingCommon");
		typedef PsimagLite::IoSelector::Out::Serializer SerializerType;
		if (counter == 0) io.createGroup(prefix);

		io.write(counter + 1,
		         prefix + "/Size",
		         (counter == 0) ? SerializerType::NO_OVERWRITE :
		                          SerializerType::ALLOW_OVERWRITE);

		prefix += ("/" + ttos(counter));

		io.createGroup(prefix);
		return prefix;
	}

	void multiSitePush(DmrgSerializerType const* ds) const
	{
		commonTargeting_.aoe().multiSitePush(ds);
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

