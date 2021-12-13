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
	typedef typename PsimagLite::Vector<RealType>::Type VectorRealType;
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
		PsimagLite::String algo = "Krylov";
		VectorRealType chebyTransform;
		SizeType disposition = 0;
		extractParamsFromName(tau, timeSteps, advanceEach, algo, chebyTransform, disposition, name);

		SizeType firstIndex = aux_.pIndexOutput();
		if (firstIndex >= aux_.pVectors().origPvectors())
			err("Cannot set an NGST to a temporary vector\n");

		AuxiliaryType* auxPtr = const_cast<AuxiliaryType*>(&aux_);

		auxPtr->pVectors().initTimeVectors(timeSteps, tau, algo, chebyTransform);

		OneTimeEvolutionType* oneTimeEvolution = aux_.timeEvolve().findThisEvolution(firstIndex);

		if (!oneTimeEvolution) {
			oneTimeEvolution = new OneTimeEvolutionType(firstIndex,
			                                            *srcVwo,
			                                            disposition,
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

private:

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
	                                  PsimagLite::String& algo,
	                                  VectorRealType& chebyTransform,
	                                  SizeType& disposition,
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
		case 5:
			if (tokens[4] != "~")
				disposition = PsimagLite::atoi(tokens[4]);
			// fall through
		case 4:
			if (tokens[3] != "~")
				algo = getChebyIfNeeded(chebyTransform, tokens[3]);
			// fall through
		case 3:
			if (tokens[2] != "~")
				advanceEach = PsimagLite::atoi(tokens[2]);
			// fall through
		case 2:
			if (tokens[1] != "~")
				timeSteps = PsimagLite::atoi(tokens[1]);
			// fall through
		case 1:
			if (tokens[0] != "~")
				tau = PsimagLite::atof(tokens[0]);
			break;
		case 0:
			break;
		default:
			err("Up to three tokens can be given\n");
			break;
		}
	}

	static PsimagLite::String getChebyIfNeeded(VectorRealType& v, PsimagLite::String algo)
	{
		static const PsimagLite::String cheby = "Chebyshev";
		if (algo.substr(0, cheby.length()) != cheby) return algo;

		VectorStringType tokens;
		PsimagLite::split(tokens, algo, ":");
		if (tokens.size() != 3 || tokens[0] != cheby)
			err("Expected Chebyshev:a:b\n");

		v.resize(2);
		v[0] = PsimagLite::atof(tokens[1]);
		v[1] = PsimagLite::atof(tokens[2]);

		return algo;
	}

	const AuxiliaryType& aux_;
	PsimagLite::ProgressIndicator progress_;
};
}
#endif // NONLOCALFORTARGETINGEXPRESSION_H
