#ifndef IMPURITYSOLVER_H
#define IMPURITYSOLVER_H

#include "PsiBase64.h"
#include "InputNg.h"
#include "LanczosSolver.h"
#include "Vector.h"
#include "ParamsDmftSolver.h"
#include "DmrgRunner.h"

namespace Dmft {

template<typename ParamsDmftSolverType>
class ImpuritySolver {

public:

	typedef typename ParamsDmftSolverType::ComplexOrRealType ComplexOrRealType;
	typedef typename PsimagLite::Real<ComplexOrRealType>::Type RealType;
	typedef typename PsimagLite::Vector<RealType>::Type VectorRealType;
	typedef DmrgRunner<ComplexOrRealType> DmrgRunnerType;
	typedef typename DmrgRunnerType::InputNgType InputNgType;

	ImpuritySolver(const ParamsDmftSolverType& params)
	    : params_(params), runner_(params_.precision)
	{}

	// bathParams[0-nBath-1] ==> V ==> hoppings impurity --> bath
	// bathParams[nBath-...] ==> energies on each bath site
	void solve(const VectorRealType& bathParams)
	{

		PsimagLite::String data;
		InputNgType::Writeable::readFile(data, params_.gsTemplate);
		PsimagLite::String data2 = modifyBathParams(data, bathParams);
		PsimagLite::String sOptions = "";

		runner_.doOneRun(data2, sOptions, params_.echoInput);
	}

	ComplexOrRealType gimp(SizeType i)
	{
		throw PsimagLite::RuntimeError("gimp not ready yet\n");
	}

private:

	static PsimagLite::String modifyBathParams(PsimagLite::String data,
	                                           const VectorRealType& bathParams)
	{
		const SizeType nBath = int(bathParams.size() / 2);
		static const PsimagLite::String label = "dir0:Connectors=";
		size_t pos1 = data.find(label);
		if (pos1 == PsimagLite::String::npos)
			err("modifyBathParams(): cannot find " + label + "\n");

		PsimagLite::String connectors = findBathParams(0, nBath, bathParams);
		PsimagLite::String buffer = data.substr(0, pos1) + label + "[" + connectors + "]";

		size_t pos2 = data.find(";", pos1);
		if (pos2 < pos1 || pos2 == PsimagLite::String::npos)
			err("modifyBathParams(): Internal error\n");

		static const PsimagLite::String label2 = "potentialV=";
		size_t pos3 = data.find(label2);
		if (pos3 < pos2 + 1 || pos3 == PsimagLite::String::npos)
			err("modifyBathParams(): cannot find " + label2 + "\n");

		SizeType len = pos3 - pos2 - 1;
		if (bathParams.size() % 2)
			err("BathParams not even\n");

		PsimagLite::String potentialV = findBathParams(nBath, 2*nBath, bathParams);
		buffer += data.substr(pos2, len) + label2 + "[0, " + potentialV +
		        ", 0, " + potentialV + "]";

		size_t pos4 = data.find(";", pos3);
		if (pos4 < pos3 + 1 || pos4 == PsimagLite::String::npos)
			err("modifyBathParams(): internal error(2)\n");

		const SizeType len2 = data.length() - pos4 - 1;
		buffer += data.substr(pos4, len2);
		return buffer;
	}

	static PsimagLite::String findBathParams(SizeType start,
	                                         SizeType end,
	                                         const VectorRealType& bathParams)
	{
		PsimagLite::String buffer = ttos(bathParams[start]);
		for (SizeType i = start + 1; i < end; ++i)
			buffer += "," + ttos(bathParams[i]);

		return buffer;
	}

	const ParamsDmftSolverType& params_;
	DmrgRunnerType runner_;
};
}
#endif // IMPURITYSOLVER_H
