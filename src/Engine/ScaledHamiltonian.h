#ifndef SCALEDHAMILTONIAN_H
#define SCALEDHAMILTONIAN_H
#include "Vector.h"
#include "ProgressIndicator.h"

namespace Dmrg {

// H' = c*H + d
// c = oneovera = (2.0-epsilon)/Wstar
// d = -oneovera*b = (2.0-epsilon)*(E0_+Wstar*0.5)/Wstar
// FIXME: TODO: Constructor should only take c and d not tstruct and E0
template<typename MatrixLanczosType, typename TargetParamsType>
class ScaledHamiltonian {

public:

	typedef typename TargetParamsType::RealType RealType;
	typedef typename MatrixLanczosType::ComplexOrRealType ComplexOrRealType;
	typedef typename PsimagLite::Vector<ComplexOrRealType>::Type VectorType;

	ScaledHamiltonian(const MatrixLanczosType& mat,
	                  const TargetParamsType& tstStruct,
	                  const RealType& E0)
	    : matx_(mat),
	      tstStruct_(tstStruct),
	      E0_(E0),
	      c_((2.0 - tstStruct_.chebyEpsilon())/tstStruct_.chebyWstar()),
	      d_(-c_*(E0_ + tstStruct_.chebyWstar()*0.5))
	{
		PsimagLite::ProgressIndicator progress("InternalMatrix");
		PsimagLite::OstringStream msg;
		msg<<"H'="<<c_<<"*H "<<d_;
		progress.printline(msg, std::cout);
	}

	SizeType rows() const { return matx_.rows(); }

	void matrixVectorProduct (VectorType &x,const VectorType &y) const
	{
		VectorType tmp(y.size());
		matx_.matrixVectorProduct(tmp, y);
		x += c_*tmp;
		x += d_*y;
	}

private:

	const MatrixLanczosType& matx_;
	const TargetParamsType& tstStruct_;
	const RealType& E0_;
	RealType c_;
	RealType d_;
}; // class ScaledHamiltonian
}

#endif // SCALEDHAMILTONIAN_H
