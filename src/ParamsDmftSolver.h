#ifndef PARAMSDMFTSOLVER_H
#define PARAMSDMFTSOLVER_H
#include "Vector.h"
#include "InputNg.h"

namespace Dmft {

template<typename ComplexOrRealType>
struct ParamsDmftSolver {

	ParamsDmftSolver()
	typedef typename PsimagLite::Real<ComplexOrRealType>::Type RealType;

	RealType ficticiousBeta = 0.01;
	SizeType nMatsubara = 100;
	SizeType numberOfKpoints = 50;
	RealType mu = 0.0;
};
}
#endif // PARAMSDMFTSOLVER_H
