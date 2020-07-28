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

	RealType operator()(const VectorRealType& args) const
	{
		RealType sum = 0.0;
		const SizeType totalMatsubaras = gammaG_.totalMatsubaras();
		for (SizeType i = 0; i < totalMatsubaras; ++i) {
			const ComplexOrRealType iwn(0, gammaG_.omega(i));
			const ComplexOrRealType val = anderson(args, iwn) - gammaG_(i);
			sum += PsimagLite::real(val*PsimagLite::conj(val));
		}

		return sum;
	}

	void df(VectorRealType& dest, const VectorRealType& src) const
	{
		for (SizeType j = 0; j < 2*nBath_; ++j) {
			RealType sum = 0.0;
			const SizeType totalMatsubaras = gammaG_.totalMatsubaras();
			for (SizeType i = 0; i < totalMatsubaras; ++i) {
				const ComplexOrRealType iwn(0, gammaG_.omega(i));
				const ComplexOrRealType val = anderson(src, iwn) - gammaG_(i);

				const ComplexOrRealType valPrime = andersonPrime(src, iwn, j);

				sum += PsimagLite::real(val*PsimagLite::conj(valPrime) +
				                        valPrime*PsimagLite::conj(val));
			}

			dest[j] = sum;
		}
	}

private:

	ComplexOrRealType anderson(const VectorRealType& args, ComplexOrRealType iwn) const
	{
		assert(args.size() == 2*nBath_);
		ComplexOrRealType sum = 0.0;
		for (SizeType i = 0; i < nBath_; ++i) {
			const RealType valpha = args[i];
			const RealType epsilon = args[i + nBath_];
			sum += valpha*valpha/(iwn - epsilon);
		}

		return sum;
	}

	ComplexOrRealType andersonPrime(const VectorRealType& args,
	                                ComplexOrRealType iwn,
	                                SizeType jnd) const
	{
		assert(args.size() == 2*nBath_);
		assert(jind < 2*nBath_);
		const RealType valpha = args[jnd];
		const RealType epsilon = args[jnd + nBath_];
		return (jnd < nBath_) ? 2.0*valpha/(iwn - epsilon) :
		                        -squareOf(valpha/(iwn - epsilon));
	}

	static ComplexOrRealType squareOf(ComplexOrRealType x) { return x*x; }

	SizeType nBath_;
	const FunctionOfFrequencyType& gammaG_;
};

}
#endif // ANDERSONFUNCTION_H
