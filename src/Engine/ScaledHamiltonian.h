#ifndef SCALEDHAMILTONIAN_H
#define SCALEDHAMILTONIAN_H
#include "Vector.h"
#include "ProgressIndicator.h"
#include "ProgramGlobals.h"

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
	typedef PsimagLite::Matrix<ComplexOrRealType> MatrixType;
	typedef typename PsimagLite::Vector<RealType>::Type VectorRealType;

	ScaledHamiltonian(const MatrixLanczosType& mat,
	                  const TargetParamsType& tstStruct,
	                  const RealType& E0,
	                  ProgramGlobals::VerboseEnum verbose)
	    : matx_(mat),
	      tstStruct_(tstStruct),
	      E0_(E0)
	{
		if (tstStruct_.chebyTransform().size() != 2)
			err("ChebyshevTransform must be a vector of two real entries\n");
		c_ = tstStruct_.chebyTransform()[0];
		d_ = tstStruct_.chebyTransform()[1];

		SizeType n = matx_.rows();
		PsimagLite::ProgressIndicator progress("InternalMatrix");
		PsimagLite::OstringStream msg;
		msg<<"H'="<<c_<<"*H "<<d_<<"      (rank="<<n<<")";
		progress.printline(msg, std::cout);

		if (verbose == ProgramGlobals::VerboseEnum::NO) return;

		MatrixType dense;
		VectorRealType eigs(n);
		matx_.fullDiag(eigs, dense);
		PsimagLite::OstringStream msg2;
		msg2<<"eigs[0]="<<eigs[0]<<" eigs["<<(n-1)<<"]="<<eigs[n - 1];
		progress.printline(msg2, std::cout);
	}

	SizeType rows() const { return matx_.rows(); }

	void matrixVectorProduct (VectorType &x,const VectorType &y) const
	{
		VectorType tmp(y.size());
		matx_.matrixVectorProduct(tmp, y);
		x += c_*tmp;
		x += d_*y;
	}

	void scale (VectorRealType& eigs) const
	{
		for (SizeType i=0; i<eigs.size();++i) {
			eigs[i] = c_*eigs[i]+d_;
		}
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
