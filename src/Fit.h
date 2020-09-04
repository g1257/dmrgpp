#ifndef FIT_H
#define FIT_H
#include "Minimizer.h"
#include "MinParams.h"
#include "AndersonFunction.h"
#include "MersenneTwister.h"

namespace Dmft {

template<typename ComplexOrRealType>
class Fit {

public:

	typedef typename PsimagLite::Real<ComplexOrRealType>::Type RealType;
	typedef typename PsimagLite::Vector<RealType>::Type VectorRealType;
	typedef MinParams<RealType> MinParamsType;
	typedef AndersonFunction<ComplexOrRealType> AndersonFunctionType;
	typedef typename AndersonFunctionType::FunctionOfFrequencyType FunctionOfFrequencyType;
	typedef PsimagLite::MersenneTwister RngType;

	Fit(SizeType nBath, const MinParamsType& minParams)
	    : nBath_(nBath), minParams_(minParams), results_(2*nBath), rng_(1234)
	{}

	// Compute the optimized bath parameters and store them in vector gammaG
	// See AndersonFunction.h documentation for the fit function, and
	// for the order of storage of bath parameters
	void fit(const FunctionOfFrequencyType& gammaG)
	{
		for (SizeType i = 0; i < results_.size(); ++i)
			results_[i] = 5.0*rng_();

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
	RngType rng_;
};
}
#endif // FIT_H
