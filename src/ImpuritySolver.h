#ifndef IMPURITYSOLVER_H
#define IMPURITYSOLVER_H
#include "Vector.h"

namespace Dmft {

template<typename ComplexOrRealType>
class ImpuritySolver {

public:

	typedef typename PsimagLite::Real<ComplexOrRealType>::Type RealType;
	typedef typename PsimagLite::Vector<RealType>::Type VectorRealType;

	void solve(const VectorRealType& bathParams) {}

	ComplexOrRealType gimp(SizeType i)
	{
		return 0.0;
	}
};
}
#endif // IMPURITYSOLVER_H
