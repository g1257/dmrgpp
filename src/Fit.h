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
		MinParams(RealType d, RealType d2, RealType t, SizeType m)
		    : delta(d), delta2(d2), tolerance(t), maxIter(m)
		{}

		RealType delta;
		RealType delta2;
		RealType tolerance;
		SizeType maxIter;
	};

	class AndersonFunction {

	};

	Fit(SizeType nBath, MinParams& minParams) : results_(2*nBath), minParams_(minParams)
	{}

	void fit()
	{
		AndersonFunction f;
		PsimagLite::Minimizer<RealType, AndersonFunction> min(f, minParams_.maxIter);
		int iter = min.conjugateGradient(results_,
		                                 minParams_.delta,
		                                 minParams_.delta2,
		                                 minParams_.tolerance);
		if (iter < 0)
            std::cerr<<"No minimum found\n";
	}

private:

	const MinParams& minParams_;
	VectorRealType results_;
};
}
#endif // FIT_H
