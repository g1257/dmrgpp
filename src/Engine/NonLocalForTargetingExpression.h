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
	typedef typename TargetParamsType::RealType RealType;
	typedef typename PsimagLite::Vector<VectorWithOffsetType>::Type VectorVectorWithOffsetType;
	typedef typename AuxiliaryType::InputValidatorType InputValidatorType;
	typedef PsimagLite::OneOperatorSpec OneOperatorSpecType;
	typedef typename OneOperatorSpecType::SiteSplit SiteSplitType;
	typedef typename ModelType::VectorSizeType VectorSizeType;
	typedef PsimagLite::Vector<PsimagLite::String>::Type VectorStringType;
	typedef std::pair<PsimagLite::String, SizeType> PairStringSizeType;

	class OneTimeEvolution {

	public:

		OneTimeEvolution(SizeType firstIndex,
		                 const VectorWithOffsetType& src,
		                 SizeType timeSteps,
		                 const AuxiliaryType& aux)
		    : indices(timeSteps)
		{
			indices[0] = firstIndex;
			AuxiliaryType* auxPtr = const_cast<AuxiliaryType*>(&aux);
			for (SizeType i = 1; i < timeSteps; ++i) {
				auto lambda = [this, i](SizeType ind) {
					this->indices[i] = ind;
					return "|P" + ttos(ind) + ">";
				};

				auxPtr->pVectors().createNew(src, lambda);
			}
		}

		VectorSizeType indices;
	};

	typedef typename PsimagLite::Vector<OneTimeEvolution*>::Type VectorOneTimeEvolutionType;

	NonLocalForTargetingExpression(const AuxiliaryType& aux)
	    : aux_(aux)
	{}

	~NonLocalForTargetingExpression()
	{
		const SizeType n = vEvolutions_.size();
		for (SizeType i = 0; i < n; ++i) {
			delete vEvolutions_[i];
			vEvolutions_[i] = nullptr;
		}
	}

	bool timeEvolve(const SiteSplitType& siteSplit,
	                PsimagLite::String srcKet,
	                SizeType site)
	{
		if (siteSplit.hasSiteString)
			err("Global operators cannot have a site\n");
		const VectorWithOffsetType& srcVwo = aux_.pVectors().getCurrentVectorConst(srcKet);
		if (srcVwo.size() == 0) return false;

		static const bool allOperatorsApplied = true;

		VectorSizeType block1(1, site);
		static const bool isLastCall = true;

		SizeType timeSteps = 3; // Fixme read from string TODO FIXME
		RealType tau = 0.1; // Fixme read from string TODO FIXME
		SizeType advanceEach = aux_.pVectors().aoe().model().superGeometry().numberOfSites() - 2;

		SizeType firstIndex = aux_.pIndexOutput();
		if (firstIndex >= aux_.pVectors().origPvectors())
			err("Cannot set an NGST to a temporary vector\n");

		AuxiliaryType* auxPtr = const_cast<AuxiliaryType*>(&aux_);

		if (vEvolutions_.size() == 0)
			auxPtr->pVectors().initTimeVectors(timeSteps, tau, advanceEach);

		OneTimeEvolution* oneTimeEvolution = findThisEvolution(firstIndex);

		if (!oneTimeEvolution) {
			oneTimeEvolution = new OneTimeEvolution(firstIndex, srcVwo, timeSteps, aux_);
			vEvolutions_.push_back(oneTimeEvolution);
		}

		auxPtr->pVectors().aoeNonConst().calcTimeVectors(oneTimeEvolution->indices,
		                                                 aux_.Eg(),
		                                                 srcVwo,
		                                                 aux_.direction(),
		                                                 allOperatorsApplied,
		                                                 false, // don't wft or advance indices[0]
		                                                 block1,
		                                                 isLastCall);
		return true;
	}

private:

	OneTimeEvolution* findThisEvolution(SizeType firstIndex) const
	{
		const SizeType n = vEvolutions_.size();
		for (SizeType i = 0; i < n; ++i)
			if (vEvolutions_[i]->indices[0] == firstIndex) return vEvolutions_[i];

		return nullptr;
	}

	const AuxiliaryType& aux_;
	static VectorOneTimeEvolutionType vEvolutions_;
};

template<typename T>
typename NonLocalForTargetingExpression<T>::VectorOneTimeEvolutionType
NonLocalForTargetingExpression<T>::vEvolutions_;
}
#endif // NONLOCALFORTARGETINGEXPRESSION_H
