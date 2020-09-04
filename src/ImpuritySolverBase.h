#ifndef IMPURITYSOLVER_BASE_H
#define IMPURITYSOLVER_BASE_H

#include "Vector.h"
#include "PsimagLite.h"

namespace Dmft {

template<typename ParamsDmftSolverType>
class ImpuritySolverBase {

public:

	typedef typename ParamsDmftSolverType::ComplexOrRealType ComplexOrRealType;
	typedef typename PsimagLite::Real<ComplexOrRealType>::Type RealType;
	typedef std::complex<RealType> ComplexType;
	typedef typename PsimagLite::Vector<RealType>::Type VectorRealType;
	typedef typename PsimagLite::Vector<ComplexType>::Type VectorComplexType;
	typedef PsimagLite::PsiApp ApplicationType;

	// bathParams[0-nBath-1] ==> V ==> hoppings impurity --> bath
	// bathParams[nBath-...] ==> energies on each bath site
	virtual void solve(const VectorRealType& bathParams) = 0;

	virtual const VectorComplexType& gimp() const  = 0;
};
}
#endif // IMPURITYSOLVER_BASE_H
