#ifndef PARAMSDMFTSOLVER_H
#define PARAMSDMFTSOLVER_H
#include "Vector.h"
#include "InputNg.h"
#include "MinParams.h"

namespace Dmft {

template<typename ComplexOrRealType_, typename InputNgType>
struct ParamsDmftSolver {

	typedef ComplexOrRealType_ ComplexOrRealType;
	typedef typename PsimagLite::Real<ComplexOrRealType>::Type RealType;
	typedef MinParams<RealType> MinParamsType;

	ParamsDmftSolver(typename InputNgType::Readable& io)
	    : echoInput(false)
	{
		io.readline(ficticiousBeta, "FicticiousBeta=");
		io.readline(mu, "ChemicalPotential=");
		io.readline(nMatsubaras, "Matsubaras=");
		io.readline(latticeGf, "LatticeGf=");
		io.readline(nBath, "NumberOfBathPoints=");
		io.readline(dmftIter, "DmftNumberOfIterations=");
		io.readline(dmftError, "DmftTolerance=");
		io.readline(gsTemplate, "DmrgGsTemplate=");
		io.readline(omegaTemplate, "DmrgOmegaTemplate=");

		try {
			io.readline(precision, "Precision=");
		} catch (std::exception&) {}

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

	bool echoInput;
	RealType ficticiousBeta;
	RealType mu;
	RealType dmftError;
	SizeType nMatsubaras;
	PsimagLite::String latticeGf;
	SizeType nBath;
	SizeType dmftIter;
	SizeType precision;
	PsimagLite::String gsTemplate;
	PsimagLite::String omegaTemplate;
	MinParamsType minParams;
};
}
#endif // PARAMSDMFTSOLVER_H
