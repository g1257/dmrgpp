#ifndef DENSITYOFSTATES_H
#define DENSITYOFSTATES_H
#include "PsimagLite.h"

namespace Dmft {

template<typename ComplexOrRealType>
class DensityOfStates {

public:

	typedef typename PsimagLite::Real<ComplexOrRealType>::Type RealType;

	DensityOfStates(PsimagLite::String option, RealType wOverTwo)
	    : wOverTwo_(wOverTwo)
	{
		if (option != "semicircular")
			err("DensityOfStates " + option +
			    " not yet supported; only semicircular supported.\n");

	}

	RealType lowerBound() const { return -wOverTwo_; }

	RealType upperBound() const { return wOverTwo_; }

	RealType operator()(RealType e) const
	{
		const RealType wOverTwoSquared = wOverTwo_*wOverTwo_;
		return 2.0*sqrt(wOverTwoSquared - e*e)/(wOverTwoSquared*M_PI);
	}

private:

	RealType wOverTwo_;
};

}
#endif // DENSITYOFSTATES_H
