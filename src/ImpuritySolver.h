#ifndef IMPURITYSOLVER_H
#define IMPURITYSOLVER_H

#include "PsiBase64.h"
#include "InputNg.h"
#include "LanczosSolver.h"
#include "Vector.h"
#include "ParamsDmftSolver.h"
#include "../../dmrgpp/src/Engine/DmrgRunner.h"
#include "../../dmrgpp/src/Engine/ManyOmegas.h"
#include "../../dmrgpp/src/Engine/ProcOmegas.h"
#include "PsimagLite.h"
#include "Matsubaras.h"

namespace Dmft {

template<typename ParamsDmftSolverType>
class ImpuritySolver {

public:

	typedef typename ParamsDmftSolverType::ComplexOrRealType ComplexOrRealType;
	typedef typename PsimagLite::Real<ComplexOrRealType>::Type RealType;
	typedef std::complex<RealType> ComplexType;
	typedef typename PsimagLite::Vector<RealType>::Type VectorRealType;
	typedef typename PsimagLite::Vector<ComplexType>::Type VectorComplexType;
	typedef Dmrg::DmrgRunner<RealType> DmrgRunnerType;
	typedef typename DmrgRunnerType::InputNgType InputNgType;
	typedef PsimagLite::PsiApp ApplicationType;
	typedef Matsubaras<RealType> MatsubarasType;
	typedef Dmrg::ManyOmegas<RealType, MatsubarasType> ManyOmegasType;
	typedef Dmrg::ProcOmegas<RealType, MatsubarasType> ProcOmegasType;

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
		PsimagLite::String insitu = "<gs|nup|gs>";

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

		const PsimagLite::String rootIname = "input";
		const PsimagLite::String rootOname = "OUTPUT";
		const bool skipFourier = true;

		Dmrg::InputCheck inputCheck;
		typename InputNgType::Writeable ioW(inputCheck, data4);
		typename InputNgType::Readable io(ioW);
		ProcOmegasType procOmegas(io,
		                          params_.precision,
		                          skipFourier,
		                          rootIname,
		                          rootOname,
		                          matsubaras);

		procOmegas.run();

		readGimp(rootOname, matsubaras);
	}

	ComplexOrRealType gimp(SizeType i)
	{
		assert(i < gimp_.size());
		return gimp_[i];
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

	void readGimp(PsimagLite::String filename, const MatsubarasType& matsubaras)
	{
		std::ifstream fin(filename);
		if (!fin || !fin.good() || fin.bad())
			err("readGimp: Could not open " + filename + "\n");

		gimp_.resize(matsubaras.total());
		SizeType ind = 0;
		while (!fin.eof()) {
			RealType val = 0;
			fin>>val;
			SizeType n = 0;
			fin>>n;
			SizeType site = 0;
			fin>>site;
			if (site != 0)
				err("readGimp: Expecting site 0, but found " + ttos(site) + " instead\n");

			RealType val1 = 0;
			fin>>val1;

			RealType val2 = 0;
			fin>>val2;

			if (ind >= gimp_.size())
				break;

			gimp_[ind++] = ComplexType(val1, val2);

			if (n == 1) continue;

			--n;
			for (SizeType i = 0; i < 3*n; ++i)
				fin>>val1;
		}

		if (ind < gimp_.size())
			err("readGimp: Not all values computed\n");
	}

	const ParamsDmftSolverType& params_;
	DmrgRunnerType runner_;
	VectorComplexType gimp_;
};
}
#endif // IMPURITYSOLVER_H
