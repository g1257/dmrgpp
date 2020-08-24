#ifndef FUNCTIONOFFREQUENCY_H
#define FUNCTIONOFFREQUENCY_H
#include "Vector.h"
#include <cassert>
#include "Matsubaras.h"

namespace Dmft {

template<typename ComplexOrRealType>
class FunctionOfFrequency {

public:

	typedef typename PsimagLite::Real<ComplexOrRealType>::Type RealType;
	typedef typename PsimagLite::Vector<ComplexOrRealType>::Type VectorType;
	typedef Matsubaras<RealType> MatsubarasType;

	FunctionOfFrequency(RealType fictBeta, SizeType nMatsubara)
	    : matsubaras_(fictBeta, nMatsubara),
	      data_(2*nMatsubara)
	{}

	// Total number of Matsubaras used (includes negatives and positives)
	SizeType totalMatsubaras() const { return matsubaras_.total(); }

	const RealType& fictitiousBeta() const { return matsubaras_.fictitiousBeta(); }

	// Matsubara number i, starts at 0, and the 0th is the most negative.
	const RealType& omega(SizeType i) const
	{
		return matsubaras_.omega(i);
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

	friend std::ostream& operator<<(std::ostream& os, const FunctionOfFrequency& f)
	{
		const SizeType n = f.data_.size();
		os<<n<<"\n";
		for (SizeType i = 0; i < n; ++i)
			os<<f.matsubaras_.omega(i)<<" "<<f.data_[i]<<"\n";
		return os;
	}

private:

	MatsubarasType matsubaras_;
	VectorType data_;           // value of this function at each point or omega
};
}
#endif // FUNCTIONOFFREQUENCY_H
