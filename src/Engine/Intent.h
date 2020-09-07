#ifndef INTENT_H
#define INTENT_H
#include "InputNg.h"
#include "AnsiColors.h"

namespace Dmrg {

template<typename ModelType>
class Intent {

	typedef PsimagLite::Vector<SizeType>::Type VectorSizeType;
	typedef typename ModelType::InputValidatorType InputValidatorType;
	typedef typename ModelType::ParametersType DmrgSolverParamsType;

	enum class IntentEnum {NONE, UNKNOWN, GS, GIMP_MATSUBARA, ARPES0, ARPES1, NEUTRONS_SZSZ};

public:

	Intent(const ModelType& model) : model_(model) {}

	void check() const
	{
		const IntentEnum intent = getIntent();

		switch (intent) {
		case IntentEnum::ARPES0:
			spectroscopy();
			break;
		case IntentEnum::ARPES1:
			spectroscopy();
			break;
		case IntentEnum::NEUTRONS_SZSZ:
			spectroscopy();
			dynTypeShouldBe(0);
			operatorShouldBe("sz");
			break;
		case IntentEnum::GS:
		case IntentEnum::GIMP_MATSUBARA:
			break;
		case IntentEnum::NONE:
			saySomethingAbout("No intent found");
			break;
		default:
			saySomethingAbout("Could not understand your intent");
			break;
		}
	}

private:

	IntentEnum getIntent() const
	{
		PsimagLite::String intent("");
		try {
			model_.ioIn().readline(intent, "Intent");
		} catch (std::exception&) {
			return IntentEnum::NONE;
		}

		if (intent == "GroundState") return IntentEnum::GS;
		if (intent == "GimpMatsubara") return IntentEnum::GIMP_MATSUBARA;
		if (intent == "ARPES0") return IntentEnum::ARPES0;
		if (intent == "ARPES1") return IntentEnum::ARPES1;
		if (intent == "NeutronsSzSz") return IntentEnum::NEUTRONS_SZSZ;
		return IntentEnum::UNKNOWN;
	}

	void spectroscopy() const
	{
		if (!hasInSolverOptions("restart"))
			saySomethingAbout(PsimagLite::String("restart") + "in SolverOptions");

		if (!hasInSolverOptions("CorrectionVectorTargeting") &&
		        !hasInSolverOptions("TargetingChebyshev"))
			saySomethingAbout(PsimagLite::String("CorrectionVectorTargeting or ") +
			                  "TargetingChebyshev in SolverOptions");

		if (!hasInSolverOptions("minimizeDisk"))
			suggest("minimizeDisk in SolverOptions");

		PsimagLite::String cvft;
		model_.ioIn().readline(cvft, "CorrectionVectorFreqType");

		if (cvft != "Real")
			saySomethingAbout("CorrectionVectorFreqType should be equal to Real");

		VectorSizeType tspsites;
		VectorSizeType tsploops;
		model_.ioIn().read(tspsites, "TSPSites");
		model_.ioIn().read(tsploops, "TSPLoops");
		if (tspsites.size() != 1 || tsploops.size() != 1) {
			saySomethingAbout("TSPSites and TSPLoops should have one site");
			return;
		}

		if (model_.params().insitu == "")
			saySomethingAbout("No insitu measurements found");

		assert(tspsites.size() > 0);
		const SizeType site = tspsites[0];
		examineSite(site);

		const SizeType finiteLoops = model_.params().finiteLoop.size();

		if (finiteLoops < tsploops.size() + 2)
			saySomethingAbout("Expected finite loops equal or bigger than tsploops + 2");
	}

	void examineSite(SizeType site) const
	{
		const SizeType l = model_.superGeometry().numberOfSites();
		const SizeType leg = (l & 1) ? 0 : getLeg();
		const SizeType lOverTwo = l/2;

		if (leg == 1) {
			if (site != lOverTwo)
				saySomethingAbout("Expected: TSPSites 1" + ttos(lOverTwo));
			return;
		} else if (leg == 2) {
			if (site != lOverTwo && site + 2 != lOverTwo)
				saySomethingAbout("Expected: TSPSites 1 L/2 or L/2 - 2");
			return;
		} else if (leg == 4) {
			if (site != lOverTwo && site != lOverTwo + 2)
				saySomethingAbout("Expected: TSPSites 1 L/2 or L/2 + 2");
			return;
		} else if (leg == 0) {
			saySomethingAbout("Couldn't check TSPSites value");
			return;
		} else {
			saySomethingAbout("Geometry may not support center site approximation");
		}
	}

	bool hasInSolverOptions(PsimagLite::String what) const
	{
		return model_.params().options.isSet(what);
	}

	void saySomethingAbout(PsimagLite::String what) const
	{
		std::cerr<<"\x1b[38;5;124m";
		std::cerr<<"WARNING: "<<what<<" (given your Intent)\n";
		std::cerr<<PsimagLite::AnsiColor::reset;

		std::cout<<"\x1b[38;5;124m";
		std::cout<<"WARNING: "<<what<<" (given your Intent)\n";
		std::cout<<PsimagLite::AnsiColor::reset;
	}

	void suggest(PsimagLite::String what) const
	{
		std::cerr<<"May I suggest "<<what<<" (given your Intent)\n";
		std::cout<<"May I suggest "<<what<<" (given your Intent)\n";
	}

	SizeType getLeg() const
	{
		SizeType terms = model_.superGeometry().terms();
		PsimagLite::String name("");
		for (SizeType t = 0; t < terms; ++t) {
			if (t > 0 && name != model_.superGeometry().label(t))
				return 0;
			name = model_.superGeometry().label(t);
		}

		if (name == "chain") return 1;

		if (name != "ladder") return 0;

		SizeType leg = 0;
		try {
			model_.ioIn().readline(leg, "LadderLeg=");
		} catch (std::exception&) {
			leg = 0;
		}

		return leg;
	}

	void dynTypeShouldBe(SizeType x) const
	{
		SizeType y = 0;
		try {
			model_.ioIn().readline(y, "DynamicDmrgType=");
		} catch (std::exception&) {
			saySomethingAbout("DynamicDmrgType= not found?");
			return;
		}

		if (x == y) return;

		saySomethingAbout("DynamicDmrgType should be " + ttos(x) +
		                  ", but " + ttos(y) + " found instead");

	}

	void operatorShouldBe(PsimagLite::String x) const
	{
		PsimagLite::String y("");
		try {
			model_.ioIn().readline(y, "OperatorExpression=");
		} catch (std::exception&) {
			saySomethingAbout("OperatorExpression= not found");
			return;
		}

		if (x == y) return;

		saySomethingAbout("OperatorExpression should be " + x +
		                  ", but " + y + " found instead");
	}

	const ModelType& model_;
};

}
#endif // INTENT_H
