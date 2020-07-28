#ifndef FUNCTIONOFFREQUENCY_H
#define FUNCTIONOFFREQUENCY_H
#include "Vector.h"

namespace Dmft {

template<typename ComplexOrRealType>
class FunctionOfFrequency {

public:

	typedef typename PsimagLite::Real<ComplexOrRealType>::Type RealType;
	typedef typename PsimagLite::Vector<RealType>::Type VectorRealType;
	typedef typename PsimagLite::Vector<ComplexOrRealType>::Type VectorType;

	FunctionOfFrequency(RealType& fictBeta, SizeType nMatsubara)
	    : fictBeta_(fictBeta), nMatsubara_(nMatsubara)
	{
		for (SizeType i = 0; i < 2*nMatsubara_; ++i) {
			int n = i - n;
			matsubaras_[i] = M_PI*n/fictBeta_;
		}
	}

private:

	RealType fictBeta_;
	SizeType nMatsubara_;
	VectorRealType matsubaras_;

};
}
#endif // FUNCTIONOFFREQUENCY_H
