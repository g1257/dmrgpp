#ifndef DISPERSION_H
#define DISPERSION_H
#include "Vector.h"
#include <cassert>

namespace Dmft {

template<typename ComplexOrRealType>
class Dispersion {

public:

	typedef typename PsimagLite::Real<ComplexOrRealType>::Type RealType;
	typedef typename PsimagLite::Vector<RealType>::Type VectorRealType;

	Dispersion(SizeType N) : ek_(N)
	{
		for (SizeType i = 0; i < N; ++i)
			ek_[i] = -2*cos(2*M_PI*i/N);
	}

	SizeType size() const { return ek_.size(); }

	const RealType& operator()(SizeType i) const
	{
		assert(i < ek_.size());
		return ek_[i];
	}

private:

	VectorRealType ek_;
};
}
#endif // DISPERSION_H
