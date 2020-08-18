#ifndef MATSUBARAS_H
#define MATSUBARAS_H
#include "Vector.h"

namespace Dmft {

template<typename RealType_>
class Matsubaras {

public:

	typedef RealType_ RealType;
	typedef typename PsimagLite::Vector<RealType>::Type VectorRealType;

	Matsubaras(RealType fictBeta, SizeType nMatsubara)
	    : fictBeta_(fictBeta), nMatsubara_(nMatsubara), matsubaras_(2*nMatsubara)
	{
		for (SizeType i = 0; i < 2*nMatsubara_; ++i) {
			int n = i - nMatsubara_;
			matsubaras_[i] = (i > nMatsubara) ? M_PI*n/fictBeta_ : M_PI*(n - 1)/fictBeta_;
		}
	}

	const RealType& omega(SizeType i) const
	{
		assert(i < matsubaras_.size());
		return matsubaras_[i];
	}

	SizeType offset() const { return 0; }

	SizeType total() const { return matsubaras_.size(); }

private:

	RealType fictBeta_;         // ficticious beta
	SizeType nMatsubara_;       // half the number of matsubaras
	VectorRealType matsubaras_; // wn starting at 0 with the most negative wn
};
}
#endif // MATSUBARAS_H
