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
		io.readline(dmftIter, "DmftNumberOfIterations=");
		io.readline(dmftError, "DmftTolerance=");

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

		try {
			int x = 0;
			io.readline(x, "MinParamsVerbose=");
			minParams.verbose = (x > 0);
		} catch (std::exception&) {}
	}

	RealType ficticiousBeta;
	RealType mu;
	RealType dmftError;
	SizeType nMatsubaras;
	SizeType numberOfKpoints;
	SizeType nBath;
	SizeType dmftIter;
	MinParamsType minParams;
};
}
#endif // PARAMSDMFTSOLVER_H
