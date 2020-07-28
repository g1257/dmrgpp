#ifndef MINPARAMS_H
#define MINPARAMS_H

namespace Dmft {

template<typename RealType>
struct MinParams {

	MinParams()
	    : delta(1e-3), delta2(1e-3), tolerance(1e-3), maxIter(100), verbose(false)
	{}

	MinParams(RealType d, RealType d2, RealType t, SizeType m, bool v)
	    : delta(d), delta2(d2), tolerance(t), maxIter(m), verbose(v)
	{}

	RealType delta;
	RealType delta2;
	RealType tolerance;
	SizeType maxIter;
	bool verbose;
};

}
#endif // MINPARAMS_H
