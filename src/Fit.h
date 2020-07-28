#ifndef FIT_H
#define FIT_H
#include "Minimizer.h"

namespace Dmft {

template<typename ComplexOrRealType>
class Fit {

public:

	typedef typename PsimagLite::Real<ComplexOrRealType>::Type RealType;
	typedef typename PsimagLite::Vector<RealType>::Type VectorRealType;

	struct MinParams {
		MinParams(RealType d, RealType d2, RealType t)
		    : delta(d), delta2(d2), tolerance(t)
		{}

		RealType delta;
		RealType delta2;
		RealType tolerance;
	};

	Fit(MinParams& minParams) : minParams_(minParams)
	{}

	void fit(SizeType maxIter)
	{
		AndersonFunctionType f;
		PsimagLite::Minimizer<RealType, AndersonFunctionType> min(f, maxIter);
		int iter = min.conjugateGradient(results_,
		                                 minParams_.delta,
		                                 minParams_.delta2,
		                                 minParams_.tolerance);
		if (iter<0)
            std::cerr<<"No minimum found\n";
	}

private:

	const MinParams& minParams_;
};
}
#endif // FIT_H
