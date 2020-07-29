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

	SizeType size() const { return 2*nBath_; }

	// Returns \sum_n |Anderson(Valpha, eAlpha, iwn) - GammaG(iwn)|^2
	// See the AndersonFunction below
	RealType operator()(const VectorRealType& args) const
	{
		RealType sum = 0.0;
		const SizeType totalMatsubaras = gammaG_.totalMatsubaras();
		for (SizeType i = 0; i < totalMatsubaras; ++i) {
			const ComplexOrRealType iwn(0, gammaG_.omega(i));
			const ComplexOrRealType val = anderson(args, iwn) - gammaG_(i);
			sum += PsimagLite::real(val*PsimagLite::conj(val));
		}

		return sum/totalMatsubaras;
	}

	// For each 0 <= j < 2*nBath, this function
	// returns the derivative of the function above with respect
	// to bath parameter j, evaluated at the bath parameters in src
	// and stores the result in dest[j]
	// for the order of bath parameters see AndersonFunction
	void df(VectorRealType& dest, const VectorRealType& src) const
	{
		const SizeType totalMatsubaras = gammaG_.totalMatsubaras();

		for (SizeType j = 0; j < 2*nBath_; ++j) {
			RealType sum = 0.0;
			for (SizeType i = 0; i < totalMatsubaras; ++i) {
				const ComplexOrRealType iwn(0, gammaG_.omega(i));
				const ComplexOrRealType val = anderson(src, iwn) - gammaG_(i);

				const ComplexOrRealType valPrime = andersonPrime(src, iwn, j);

				sum += PsimagLite::real(val*PsimagLite::conj(valPrime) +
				                        valPrime*PsimagLite::conj(val));
			}

			dest[j] = sum/totalMatsubaras;
		}
	}

	// Returns \sum_{0<=j<nBath} V_j^2/(iwn - epsilon_j),
	// where the V_j are stored in the first half or args,
	// and the epsilon_j are stored in the last half or args
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

	static ComplexOrRealType squareOf(ComplexOrRealType x) { return x*x; }

private:

	// For any 0 <= jnd < 2*nBath, this function returns the derivative of
	// the AndersonFunction above with respect to bath parameter jnd,
	// evaluated at the bath parameters args
	// The order in which bath parameters are stored is described
	// under AndersonFunction
	ComplexOrRealType andersonPrime(const VectorRealType& args,
	                                ComplexOrRealType iwn,
	                                SizeType jnd) const
	{
		assert(args.size() == 2*nBath_);
		assert(jnd < 2*nBath_);
		const RealType valpha = args[jnd];
		const RealType epsilon = args[jnd + nBath_];
		return (jnd < nBath_) ? 2.0*valpha/(iwn - epsilon) :
		                        squareOf(valpha/(iwn - epsilon));
	}

	SizeType nBath_; // Number of Bath sites
	const FunctionOfFrequencyType& gammaG_; // Gamma Green Function to fit
};

}
#endif // ANDERSONFUNCTION_H
