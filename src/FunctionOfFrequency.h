#ifndef FUNCTIONOFFREQUENCY_H
#define FUNCTIONOFFREQUENCY_H
#include "Vector.h"
#include <cassert>

namespace Dmft {

template<typename ComplexOrRealType>
class FunctionOfFrequency {

public:

	typedef typename PsimagLite::Real<ComplexOrRealType>::Type RealType;
	typedef typename PsimagLite::Vector<RealType>::Type VectorRealType;
	typedef typename PsimagLite::Vector<ComplexOrRealType>::Type VectorType;

	FunctionOfFrequency(RealType fictBeta, SizeType nMatsubara)
	    : fictBeta_(fictBeta),
	      nMatsubara_(nMatsubara),
	      matsubaras_(2*nMatsubara),
	      data_(2*nMatsubara)
	{
		for (SizeType i = 0; i < 2*nMatsubara_; ++i) {
			int n = i - nMatsubara_;
			matsubaras_[i] = (i > nMatsubara) ? M_PI*n/fictBeta_ : M_PI*(n - 1)/fictBeta_;
		}
	}

	// Total number of Matsubaras used (includes negatives and positives)
	SizeType totalMatsubaras() const { return matsubaras_.size(); }

	// Matsubara number i, starts at 0, and the 0th is the most negative.
	const RealType& omega(SizeType i) const
	{
		assert(i < totalMatsubaras());
		return matsubaras_[i];
	}

	// Returns the content of this function at point i
	// the wn at this point is given by omega(i) above
	// for reading value only
	const ComplexOrRealType& operator()(SizeType i) const
	{
		assert(i < totalMatsubaras());
		return data_[i];
	}

	// Same as above, but allows for modification of the value
	// for writing the value
	ComplexOrRealType& operator()(SizeType i)
	{
		assert(i < totalMatsubaras());
		return data_[i];
	}

private:

	RealType fictBeta_;         // ficticious beta
	SizeType nMatsubara_;       // half the number of matsubaras
	VectorRealType matsubaras_; // wn starting at 0 with the most negative wn
	VectorType data_;           // value of this function at each point or omega
};
}
#endif // FUNCTIONOFFREQUENCY_H
