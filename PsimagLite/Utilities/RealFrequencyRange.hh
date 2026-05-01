#ifndef REAL_FREQ_RANGE_H
#define REAL_FREQ_RANGE_H
#include "InputNg.h"
#include "Vector.h"
#include <cassert>

namespace PsimagLite {

template <typename RealType_> class RealFrequencyRange {

public:

	using RealType       = RealType_;
	using VectorRealType = typename PsimagLite::Vector<RealType>::Type;

	RealFrequencyRange(const RealType& begin,
	                   const RealType& step,
	                   SizeType        total,
	                   const RealType& delta)
	    : begin_(begin)
	    , step_(step)
	    , delta_(delta)
	    , total_(total)
	{ }

	template <typename SomeReadableType> explicit RealFrequencyRange(SomeReadableType& io)
	{
		io.readline(begin_, "OmegaBegin=");
		io.readline(step_, "OmegaStep=");
		io.readline(total_, "OmegaTotal=");
		io.readline(delta_, "OmegaDelta=");
	}

	RealType omega(SizeType i) const { return begin_ + step_ * i; }

	SizeType total() const { return total_; }

	RealType delta() const { return delta_; }

	RealType offset() const { return 0; }

private:

	RealType begin_;
	RealType step_;
	RealType delta_;
	SizeType total_;
};
}
#endif // REAL_FREQ_RANGE_H
