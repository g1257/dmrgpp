#ifndef PARAMSDMFTSOLVER_H
#define PARAMSDMFTSOLVER_H
#include "Vector.h"
#include "InputNg.h"
#include "MinParams.h"

namespace Dmft {

template<typename ComplexOrRealType, typename InputNgType>
struct ParamsDmftSolver {

	typedef typename PsimagLite::Real<ComplexOrRealType>::Type RealType;
	typedef MinParams<RealType> MinParamsType;

	ParamsDmftSolver(typename InputNgType::Readable& io)
	{
		io.readline(ficticiousBeta, "FicticiousBeta=");
		io.readline(mu, "ChemicalPotential=");
		io.readline(nMatsubaras, "Matsubaras=");
		io.readline(numberOfKpoints, "NumberOfKpoints=");
		io.readline(nBath, "NumberOfBathPoints=");

		try {
			io.readline(minParams.delta, "MinParamsDelta=");
		} catch (std::exception&) {}

		try {
			io.readline(minParams.delta2, "MinParamsDelta2=");
		} catch (std::exception&) {}

		try {
			io.readline(minParams.tolerance, "MinParamsTolerance=");
		} catch (std::exception&) {}

		try {
			io.readline(minParams.maxIter, "MinParamsMaxIter=");}
		catch (std::exception&) {}
	}

	RealType ficticiousBeta;
	RealType mu;
	SizeType nMatsubaras;
	SizeType numberOfKpoints;
	SizeType nBath;
	MinParamsType minParams;
};
}
#endif // PARAMSDMFTSOLVER_H
