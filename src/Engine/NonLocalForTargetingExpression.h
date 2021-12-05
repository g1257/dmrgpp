#ifndef NONLOCALFORTARGETINGEXPRESSION_H
#define NONLOCALFORTARGETINGEXPRESSION_H
#include "AuxForTargetingExpression.h"
#include "TargetParamsBase.h"
#include "OneOperatorSpec.h"

namespace Dmrg {

template<typename TargetingBaseType>
class NonLocalForTargetingExpression {

public:

	typedef AuxForTargetingExpression<TargetingBaseType> AuxiliaryType;
	typedef typename TargetingBaseType::VectorWithOffsetType VectorWithOffsetType;
	typedef typename TargetingBaseType::ModelType ModelType;
	typedef typename TargetingBaseType::LanczosSolverType LanczosSolverType;
	typedef typename TargetingBaseType::WaveFunctionTransfType WaveFunctionTransfType;
	typedef typename TargetingBaseType::TargetParamsType TargetParamsType;
	typedef typename PsimagLite::Vector<VectorWithOffsetType>::Type VectorVectorWithOffsetType;
	typedef typename AuxiliaryType::InputValidatorType InputValidatorType;
	typedef PsimagLite::OneOperatorSpec OneOperatorSpecType;
	typedef typename OneOperatorSpecType::SiteSplit SiteSplitType;
	typedef typename ModelType::VectorSizeType VectorSizeType;
	typedef PsimagLite::Vector<PsimagLite::String>::Type VectorStringType;

	NonLocalForTargetingExpression(const AuxiliaryType& aux)
	    : aux_(aux)
	{}

	void timeEvolve(const SiteSplitType& siteSplit,
	                PsimagLite::String srcKet,
	                PsimagLite::String mangledOp,
	                SizeType site)
	{
		if (siteSplit.hasSiteString)
			err("Global operators cannot have a site\n");
		const VectorWithOffsetType& srcVwo = aux_.getCurrentVectorConst(srcKet);

		static const bool allOperatorsApplied = true;

		VectorSizeType block1(1, site);
		static const bool isLastCall = true;

		SizeType timeSteps = 5;
		VectorSizeType indices(timeSteps);
		SizeType lastIndex = aux_.aoe().targetVectors().size();
		aoeTable(indices, lastIndex);
		AuxiliaryType* auxPtr = const_cast<AuxiliaryType*>(&aux_);
		auxPtr->aoeNonConst().calcTimeVectors(indices,
		                                      aux_.Eg(),
		                                      srcVwo,
		                                      aux_.direction(),
		                                      allOperatorsApplied,
		                                      false, // don't wft or advance indices[0]
		                                      block1,
		                                      isLastCall);
	}

private:

	void aoeTable(VectorSizeType& indices, SizeType lastIndex) // lastindex = 2
	{
		SizeType timeSteps = indices.size();
		AuxiliaryType* auxPtr = const_cast<AuxiliaryType*>(&aux_);
		PsimagLite::String root = "P" + ttos(lastIndex);
		aoeTable_.resize(lastIndex + timeSteps);
		for (SizeType i = 0; i < timeSteps; ++i) {
			aoeTable_[lastIndex + i] = root + "T" + ttos(i);
			indices[i] = auxPtr->aoeNonConst().createPvector(aux_.tempVectors()[0]);
		}
	}

	const AuxiliaryType& aux_;
	VectorStringType aoeTable_;
};
}
#endif // NONLOCALFORTARGETINGEXPRESSION_H
