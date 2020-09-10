#ifndef FIT_H
#define FIT_H
#include "Minimizer.h"
#include "MinParams.h"
#include "AndersonFunction.h"
#include "MersenneTwister.h"
#include "PsimagLite.h"

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

	struct InitResults {

		InitResults(RealType ra_,
		            RealType rb_,
		            const VectorRealType& result_,
		            bool reset_)
		    : ra(ra_), rb(rb_), result(result_), reset(reset_)
		{}

		template<typename ReadableType>
		InitResults(ReadableType& io) : ra(0), rb(0), reset(false)
		{
			try {
				io.readline(ra, "InitBathRa=");
			} catch (std::exception&) {}

			try {
				io.readline(rb, "InitBathRb=");
			} catch (std::exception&) {}

			try {
				io.read(result, "InitBathVector");
			} catch (std::exception&) {}

			try {
				int tmp = 0;
				io.readline(tmp, "InitBathReset=");
				reset = (tmp > 0);
			} catch (std::exception&) {}

			bool nonConstant = (result.size() > 0);
			bool isConstant = (ra != 0 || rb != 0);
			if (isConstant && nonConstant)
				err("InitBath: Cannot have ra or rb and also a vector of init results\n");

			if (!nonConstant && !isConstant) ra = 1;
		}

		RealType ra;
		RealType rb;
		VectorRealType result;
		bool reset;
	};

	Fit(SizeType nBath,
	    const MinParamsType& minParams,
	    const InitResults& initResults)
	    : nBath_(nBath),
	      minParams_(minParams),
	      results_(2*nBath),
	      rng_(1234),
	      initResults_(initResults)
	{
		setResults();
	}

	// Compute the optimized bath parameters and store them in vector gammaG
	// See AndersonFunction.h documentation for the fit function, and
	// for the order of storage of bath parameters
	void fit(const FunctionOfFrequencyType& gammaG)
	{
		if (initResults_.reset) setResults();

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

	void setResults()
	{
		bool nonConstant = (initResults_.result.size() > 0);
		bool isConstant = (initResults_.ra != 0 || initResults_.rb != 0);
		if (isConstant && nonConstant)
			err("InitResults: Cannot have ra or rb and also a vector of init results\n");

		const RealType constant = (isConstant) ? initResults_.ra*rng_() + initResults_.rb
		                                       : 0;

		if (nonConstant && initResults_.result.size() != results_.size())
			err(PsimagLite::String("InitResults: vector of init results has wrong size: ") +
			    "expected " + ttos(results_.size()) + ", but found " +
			    ttos(initResults_.result.size()) + "\n");

		for (SizeType i = 0; i < results_.size(); ++i)
			results_[i] = (isConstant) ? constant : initResults_.result[i];
	}

	const SizeType nBath_;                  // number of bath sites
	const MinParamsType& minParams_; // parameters for fitting algorithm
	VectorRealType results_;         // stores bath parameters
	RngType rng_;
	const InitResults& initResults_;
};
}
#endif // FIT_H
