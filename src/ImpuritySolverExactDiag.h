#ifndef IMPURITYSOLVER_EXACTD_H
#define IMPURITYSOLVER_EXACTD_H

#include "Vector.h"
#include "PsimagLite.h"
#include "ImpuritySolverBase.h"

namespace Dmft {

template<typename ParamsDmftSolverType>
class ImpuritySolverExactDiag : public ImpuritySolverBase<ParamsDmftSolverType> {

public:

	typedef typename ParamsDmftSolverType::ComplexOrRealType ComplexOrRealType;
	typedef typename PsimagLite::Real<ComplexOrRealType>::Type RealType;
	typedef std::complex<RealType> ComplexType;
	typedef typename PsimagLite::Vector<RealType>::Type VectorRealType;
	typedef typename PsimagLite::Vector<ComplexType>::Type VectorComplexType;
	typedef typename ImpuritySolverBase<ParamsDmftSolverType>::ApplicationType ApplicationType;

	ImpuritySolverExactDiag(const ParamsDmftSolverType& params, const ApplicationType& app)
	    : params_(params)
	{}

	// bathParams[0-nBath-1] ==> V ==> hoppings impurity --> bath
	// bathParams[nBath-...] ==> energies on each bath site
	void solve(const VectorRealType& bathParams)
	{
		err("Unimplemented\n");
	}

	const VectorComplexType& gimp() const { return gimp_; }

private:

	const ParamsDmftSolverType& params_;
	VectorComplexType gimp_;
};
}
#endif // IMPURITYSOLVER_EXACTD_H
