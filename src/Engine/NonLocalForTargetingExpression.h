#ifndef NONLOCALFORTARGETINGEXPRESSION_H
#define NONLOCALFORTARGETINGEXPRESSION_H
#include "AuxForTargetingExpression.h"
#include "TargetParamsBase.h"
#include "OneOperatorSpec.h"
#include "ProgressIndicator.h"

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
	typedef typename TargetingBaseType::ApplyOperatorExpressionType ApplyOperatorExpressionType;
	typedef typename PsimagLite::Vector<VectorWithOffsetType>::Type VectorVectorWithOffsetType;
	typedef typename AuxiliaryType::InputValidatorType InputValidatorType;
	typedef typename AuxiliaryType::GroupOfOneTimeEvolutionsType GroupOfOneTimeEvolutionsType;
	typedef typename GroupOfOneTimeEvolutionsType::OneTimeEvolutionType OneTimeEvolutionType;
	typedef PsimagLite::OneOperatorSpec OneOperatorSpecType;
	typedef typename OneOperatorSpecType::SiteSplit SiteSplitType;
	typedef typename ModelType::VectorSizeType VectorSizeType;
	typedef PsimagLite::Vector<PsimagLite::String>::Type VectorStringType;
	typedef std::pair<PsimagLite::String, SizeType> PairStringSizeType;
	typedef typename ApplyOperatorExpressionType::TimeVectorsBaseType TimeVectorsBaseType;

	NonLocalForTargetingExpression(const AuxiliaryType& aux)
	    : aux_(aux), progress_("NonLocalForTargetingExpression")
	{}

	bool timeEvolve(PsimagLite::String name,
	                const SiteSplitType& siteSplit,
	                PsimagLite::String srcKet,
	                SizeType site)
	{
		if (siteSplit.hasSiteString)
			err("Global operators cannot have a site\n");

		const VectorWithOffsetType* srcVwo = &aux_.pVectors().getCurrentVectorConst(srcKet);
		if (srcVwo->size() == 0) return false;

		static const bool allOperatorsApplied = true;

		VectorSizeType block1(1, site);
		static const bool isLastCall = true;

		SizeType timeSteps = 3; // Fixme read from string TODO FIXME
		RealType tau = 0.1; // Fixme read from string TODO FIXME
		SizeType advanceEach = aux_.pVectors().aoe().model().superGeometry().numberOfSites() - 2;
		extractParamsFromName(tau, timeSteps, advanceEach, name);

		SizeType firstIndex = aux_.pIndexOutput();
		if (firstIndex >= aux_.pVectors().origPvectors())
			err("Cannot set an NGST to a temporary vector\n");

		AuxiliaryType* auxPtr = const_cast<AuxiliaryType*>(&aux_);

		auxPtr->pVectors().initTimeVectors(timeSteps, tau);

		OneTimeEvolutionType* oneTimeEvolution = aux_.timeEvolve().findThisEvolution(firstIndex);

		if (!oneTimeEvolution) {
			oneTimeEvolution = new OneTimeEvolutionType(firstIndex,
			                                            *srcVwo,
			                                            timeSteps,
			                                            aux_.pVectors());
			aux_.timeEvolve().pushBack(oneTimeEvolution);
		}

		bool timeHasAdvanced = advanceInTimeOrNot(*oneTimeEvolution, advanceEach, site, tau);

		assert(oneTimeEvolution->indices().size() > 1);
		const SizeType last = oneTimeEvolution->indices().size() - 1;
		const SizeType advanceOrNot = (timeHasAdvanced) ? last : 0;
		const SizeType firstOrLast = oneTimeEvolution->indices()[advanceOrNot];
		const VectorWithOffsetType* phi = (oneTimeEvolution->time() > 0)
		        ? new VectorWithOffsetType(aux_.pVectors().aoe().targetVectors(firstOrLast))
		        : srcVwo;

		auxPtr->pVectors().aoeNonConst().calcTimeVectors(oneTimeEvolution->indices(),
		                                                 aux_.Eg(),
		                                                 *phi,
		                                                 aux_.direction(),
		                                                 allOperatorsApplied,
		                                                 false, // don't wft or advance indices[0]
		                                                 block1,
		                                                 isLastCall);

		if (oneTimeEvolution->time() > 0) {
			delete phi;
			phi = nullptr;
		}

		return true;
	}

	static bool isGlobalOperator(PsimagLite::String opName)
	{
		return isTimeEvolution(opName);
	}

private:

	static bool isTimeEvolution(PsimagLite::String op)
	{
		static const PsimagLite::String timeEvolve = "TimeEvolve";

		return (op.substr(0, timeEvolve.length()) == timeEvolve);
	}

	bool advanceInTimeOrNot(OneTimeEvolutionType& oneTimeEvolution,
	                        SizeType advanceEach,
	                        SizeType site,
	                        RealType tau)
	{
		static const bool advanceOnlyAtBorder = 1;
		const SizeType sites = aux_.pVectors().aoe().model().superGeometry().numberOfSites();
		const bool weAreAtBorder = (site == 1 || site == sites - 1);
		const bool dontAdvance = (advanceOnlyAtBorder & !weAreAtBorder);
		bool timeHasAdvanced = false;
		if (advanceEach > 0 &&
		        oneTimeEvolution.timesWithoutAdvancement() >= advanceEach
		        && !dontAdvance) {
			oneTimeEvolution.resetTimesWithoutAdvancement();
			oneTimeEvolution.advanceTime(tau);
			timeHasAdvanced = true;
		} else {
			oneTimeEvolution.incrementTimesWithoutAdvancement();
		}

		// make sure aoe.timeVectors.time_ is in sync with oneTimeEvolution's time
		TimeVectorsBaseType* ptr = const_cast<TimeVectorsBaseType*>
		        (&aux_.pVectors().aoe().timeVectors());
		ptr->setCurrentTime(oneTimeEvolution.time());

		PsimagLite::OstringStream msgg2(std::cout.precision());
		PsimagLite::OstringStream::OstringStreamType& msg2 = msgg2();
		msg2<<"Steps without advance: "<<oneTimeEvolution.timesWithoutAdvancement();
		msg2<<" site="<<site<<" currenTime="<<oneTimeEvolution.time();
		if (oneTimeEvolution.timesWithoutAdvancement() > 0)
			progress_.printline(msgg2, std::cout);
		return timeHasAdvanced;
	}

	static void extractParamsFromName(RealType& tau,
	                                  SizeType& timeSteps,
	                                  SizeType& advanceEach,
	                                  PsimagLite::String name)
	{
		//TimeEvolve{0.1,5,14}
		const PsimagLite::String tev = "TimeEvolve";
		PsimagLite::String te = name.substr(0, tev.length());
		assert(te == tev);
		name = name.substr(tev.length(), name.length() - tev.length());
		PsimagLite::String buffer;
		const SizeType m = name.length();
		for (SizeType i = 0; i < m; ++i) {
			if (name[i] == ' ' || name[i] == '{' || name[i] == '}' ||
			        name[i] == '(' || name[i] == ')') continue;
			buffer += name[i];
		}

		VectorStringType tokens;
		PsimagLite::split(tokens, buffer, ",");
		const SizeType n = tokens.size();
		switch (n) {
		case 3:
			advanceEach = PsimagLite::atoi(tokens[2]);
			// fall through
		case 2:
			timeSteps = PsimagLite::atoi(tokens[1]);
			// fall through
		case 1:
			tau = PsimagLite::atof(tokens[0]);
			break;
		case 0:
			break;
		default:
			err("Up to three tokens can be given\n");
			break;
		}
	}

	const AuxiliaryType& aux_;
	PsimagLite::ProgressIndicator progress_;
};
}
#endif // NONLOCALFORTARGETINGEXPRESSION_H
