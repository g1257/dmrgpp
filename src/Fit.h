#ifndef FIT_H
#define FIT_H
#include "Minimizer.h"
#include "MinParams.h"

namespace Dmft {

template<typename ComplexOrRealType>
class Fit {

public:

	typedef typename PsimagLite::Real<ComplexOrRealType>::Type RealType;
	typedef typename PsimagLite::Vector<RealType>::Type VectorRealType;
	typedef MinParams<RealType> MinParamsType;

	class AndersonFunction {

	};

	Fit(SizeType nBath, const MinParamsType& minParams) : minParams_(minParams), results_(2*nBath)
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

	const MinParamsType& minParams_;
	VectorRealType results_;
};
}
#endif // FIT_H
