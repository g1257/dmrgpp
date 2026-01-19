#ifndef MINPARAMS_H
#define MINPARAMS_H
#include "Vector.h"

namespace Dmft {

template <typename RealType> struct MinParams {

	enum class Method
	{
		SIMPLEX,
		CONJUGATE_GRADIENT
	};

	template <typename Readable>
	MinParams(Readable& io)
	    : delta(1e-3)
	    , delta2(1e-3)
	    , tolerance(1e-3)
	    , maxIter(100)
	    , verbose(false)
	    , method(Method::CONJUGATE_GRADIENT)
	{
		try {
			io.readline(delta, "MinParamsDelta=");
		} catch (std::exception&) { }

		try {
			io.readline(delta2, "MinParamsDelta2=");
		} catch (std::exception&) { }

		try {
			io.readline(tolerance, "MinParamsTolerance=");
		} catch (std::exception&) { }

		try {
			io.readline(maxIter, "MinParamsMaxIter=");
		} catch (std::exception&) { }

		try {
			int x = 0;
			io.readline(x, "MinParamsVerbose=");
			verbose = (x > 0);
		} catch (std::exception&) { }

		try {
			PsimagLite::String method_;
			io.readline(method_, "FitMethod=");
			method
			    = (method_ == "simplex") ? Method::SIMPLEX : Method::CONJUGATE_GRADIENT;
		} catch (std::exception&) { }
	}

	//	MinParams(RealType d,
	//	          RealType d2,
	//	          RealType t,
	//	          SizeType m,
	//	          bool v,
	//	          PsimagLite::String method_)
	//	    : delta(d), delta2(d2), tolerance(t), maxIter(m), verbose(v)
	//	{
	//		method = (method_ == "simplex") ? Method::SIMPLEX :
	//Method::CONJUGATE_GRADIENT;
	//	}

	RealType delta;
	RealType delta2;
	RealType tolerance;
	SizeType maxIter;
	bool     verbose;
	Method   method;
};

}
#endif // MINPARAMS_H
