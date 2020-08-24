#ifndef FIT_H
#define FIT_H
#include "Minimizer.h"
#include "MinParams.h"
#include "AndersonFunction.h"

namespace Dmft {

template<typename ComplexOrRealType>
class Fit {

public:

	typedef typename PsimagLite::Real<ComplexOrRealType>::Type RealType;
	typedef typename PsimagLite::Vector<RealType>::Type VectorRealType;
	typedef MinParams<RealType> MinParamsType;
	typedef AndersonFunction<ComplexOrRealType> AndersonFunctionType;
	typedef typename AndersonFunctionType::FunctionOfFrequencyType FunctionOfFrequencyType;

	Fit(SizeType nBath, const MinParamsType& minParams)
	    : nBath_(nBath), minParams_(minParams), results_(2*nBath)
	{
		for (SizeType i = 0; i < 2*nBath_; ++i) results_[i] = 0.1;
	}

	// Compute the optimized bath parameters and store them in vector gammaG
	// See AndersonFunction.h documentation for the fit function, and
	// for the order of storage of bath parameters
	void fit(const FunctionOfFrequencyType& gammaG)
	{
		AndersonFunctionType f(nBath_, gammaG);
		PsimagLite::Minimizer<RealType, AndersonFunctionType> min(f,
		                                                          minParams_.maxIter,
		                                                          minParams_.verbose);
		int iter = min.conjugateGradient(results_,
		                                 minParams_.delta,
		                                 minParams_.delta2,
		                                 minParams_.tolerance);
		if (iter < 0)
            std::cerr<<"No minimum found\n";
	}

	const VectorRealType& result() const { return results_; }

	SizeType nBath() const { return nBath_; }

private:

	SizeType nBath_;                  // number of bath sites
	const MinParamsType& minParams_; // parameters for fitting algorithm
	VectorRealType results_;         // stores bath parameters
};
}
#endif // FIT_H
