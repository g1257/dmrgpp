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

	FunctionOfFrequency(RealType& fictBeta, SizeType nMatsubara)
	    : fictBeta_(fictBeta),
	      nMatsubara_(nMatsubara),
	      matsubaras_(2*nMatsubara),
	      data_(2*nMatsubara)
	{
		for (SizeType i = 0; i < 2*nMatsubara_; ++i) {
			int n = i - nMatsubara_;
			matsubaras_[i] = M_PI*n/fictBeta_;
		}
	}

	SizeType totalMatsubaras() const { return matsubaras_.size(); }

	const RealType& omega(SizeType i) const
	{
		assert(i < totalMatsubaras());
		return matsubaras_[i];
	}

	const ComplexOrRealType& operator()(SizeType i) const
	{
		assert(i < totalMatsubaras());
		return data_[i];
	}

	ComplexOrRealType& operator()(SizeType i)
	{
		assert(i < totalMatsubaras());
		return data_[i];
	}

private:

	RealType fictBeta_;
	SizeType nMatsubara_;
	VectorRealType matsubaras_;
	VectorType data_;
};
}
#endif // FUNCTIONOFFREQUENCY_H
