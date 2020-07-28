#ifndef ANDERSONFUNCTION_H
#define ANDERSONFUNCTION_H
#include "Vector.h"
#include "FunctionOfFrequency.h"

namespace Dmft {

template<typename ComplexOrRealType>
class AndersonFunction {

public:

	typedef typename PsimagLite::Real<ComplexOrRealType>::Type RealType;
	typedef typename PsimagLite::Vector<RealType>::Type VectorRealType;
	typedef FunctionOfFrequency<ComplexOrRealType> FunctionOfFrequencyType;

	typedef RealType FieldType;

	AndersonFunction(SizeType nBath, const FunctionOfFrequencyType& gammaG)
	    : nBath_(nBath), gammaG_(gammaG)
	{}

	SizeType size() const { return nBath_; }

	RealType operator()(VectorRealType& args)
	{
		RealType sum = 0.0;
		const SizeType totalMatsubaras = gammaG_.totalMatsubaras();
		for (SizeType i = 0; i < totalMatsubaras; ++i) {
			const ComplexOrRealType iwn(0, gammaG_.omega(i));
			const ComplexOrRealType val = gammaG_(i) - anderson(args, iwn);
			sum += PsimagLite::real(val*PsimagLite::conj(val));
		}

		return sum;
	}

private:

	ComplexOrRealType anderson(VectorRealType& args, ComplexOrRealType iwn)
	{
		assert(args.size() == 2*nBath_);
		ComplexOrRealType sum = 0.0;
		for (SizeType i = 0; i < nBath_; ++i) {
			const RealType valpha = args[i];
			const RealType epsilon = args[i + nBath_];
			sum += valpha*PsimagLite::conj(valpha)/(iwn - epsilon);
		}

		return sum;
	}

	SizeType nBath_;
	const FunctionOfFrequencyType& gammaG_;
};

}
#endif // ANDERSONFUNCTION_H
