#ifndef IMPURITYSOLVER_H
#define IMPURITYSOLVER_H

#include "PsiBase64.h"
#include "InputNg.h"
#include "LanczosSolver.h"
#include "Vector.h"
#include "ParamsDmftSolver.h"
#include "../../dmrgpp/src/Engine/DmrgRunner.h"
#include "../../dmrgpp/src/Engine/ManyOmegas.h"
#include "PsimagLite.h"
#include "Matsubaras.h"

namespace Dmft {

template<typename ParamsDmftSolverType>
class ImpuritySolver {

public:

	typedef typename ParamsDmftSolverType::ComplexOrRealType ComplexOrRealType;
	typedef typename PsimagLite::Real<ComplexOrRealType>::Type RealType;
	typedef typename PsimagLite::Vector<RealType>::Type VectorRealType;
	typedef Dmrg::DmrgRunner<ComplexOrRealType> DmrgRunnerType;
	typedef typename DmrgRunnerType::InputNgType InputNgType;
	typedef PsimagLite::PsiApp ApplicationType;
	typedef Dmrg::ManyOmegas<RealType, Matsubaras<RealType> > ManyOmegasType;

	ImpuritySolver(const ParamsDmftSolverType& params, const ApplicationType& app)
	    : params_(params), runner_(params_.precision, app)
	{}

	// bathParams[0-nBath-1] ==> V ==> hoppings impurity --> bath
	// bathParams[nBath-...] ==> energies on each bath site
	void solve(const VectorRealType& bathParams)
	{

		PsimagLite::String data;
		InputNgType::Writeable::readFile(data, params_.gsTemplate);
		PsimagLite::String data2 = modifyBathParams(data, bathParams);
		PsimagLite::String insitu = "";

		runner_.doOneRun(data2, insitu, "-");

		PsimagLite::String data3;
		InputNgType::Writeable::readFile(data3, params_.omegaTemplate);
		PsimagLite::String data4 = modifyBathParams(data3, bathParams);
		PsimagLite::String insitu2 = "<gs|c'|P2>,<gs|c'|P3>";

		Matsubaras<RealType> matsubaras(params_.ficticiousBeta, params_.nMatsubaras);

		ManyOmegasType manyOmegas(data4,
		                          params_.precision,
		                          matsubaras,
		                          runner_.application());

		const bool dryrun = false;
		const PsimagLite::String rootname = "dmftDynamics";
		manyOmegas.run(dryrun, rootname, insitu2);
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

	static PsimagLite::String replaceOmega(PsimagLite::String data, RealType wn)
	{
		const PsimagLite::String omega = "$omega";
		size_t pos1 = data.find(omega, 0);
		size_t pos2 = pos1 + omega.length();
		size_t len2 = data.length() - pos1;
		return data.substr(0, pos1) + ttos(wn) + data.substr(pos2, len2);
	}

	const ParamsDmftSolverType& params_;
	DmrgRunnerType runner_;
};
}
#endif // IMPURITYSOLVER_H
