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

	enum class DmrgType {TYPE_0, TYPE_1};

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
		PsimagLite::String data2 = addBathParams(data, bathParams);
		PsimagLite::String insitu = "<gs|nup|gs>";

		runner_.doOneRun(data2, insitu, "-");

		PsimagLite::String data3;
		InputNgType::Writeable::readFile(data3, params_.omegaTemplate);
		PsimagLite::String data4 = addBathParams(data3, bathParams);

		doType(DmrgType::TYPE_0, data4);

		doType(DmrgType::TYPE_1, data4);
	}

	const VectorComplexType& gimp() const { return gimp_; }

private:

	static PsimagLite::String addBathParams(PsimagLite::String data,
	                                           const VectorRealType& bathParams)
	{
		const SizeType nBath = int(bathParams.size() / 2);
		const PsimagLite::String connectors = findBathParams(0, nBath, bathParams);
		const PsimagLite::String label = "dir0:Connectors=[" + connectors + "];\n";
		const PsimagLite::String potentialV = findBathParams(nBath, 2*nBath, bathParams);
		const PsimagLite::String label2 = "potentialV=[0, " + potentialV +
		        ", 0, " + potentialV + "];\n";

		return data + label + label2;
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

	static PsimagLite::String addTypeAndObs(DmrgType t, PsimagLite::String data)
	{
		const PsimagLite::String obsTc = (t == DmrgType::TYPE_0) ? "c" : "c'";
		const SizeType tt = (t == DmrgType::TYPE_0) ? 0 : 1;
		return data +  "integer DynamicDmrgType=" + ttos(tt) + ";\n" +
		        "OperatorExpression=\"" + obsTc + "\";\n";
	}

	void doType(DmrgType t, PsimagLite::String data)
	{
		PsimagLite::String obs = (t == DmrgType::TYPE_0) ? "c'" : "c";
		PsimagLite::String insitu2 = "<gs|" + obs + "|P2>,<gs|" + obs + "|P3>";

		PsimagLite::String data2 = addTypeAndObs(t, data);

		Matsubaras<RealType> matsubaras(params_.ficticiousBeta, params_.nMatsubaras);

		ManyOmegasType manyOmegas(data2,
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
		typename InputNgType::Writeable ioW(inputCheck, data2);
		typename InputNgType::Readable io(ioW);
		ProcOmegasType procOmegas(io,
		                          params_.precision,
		                          skipFourier,
		                          rootIname,
		                          rootOname,
		                          matsubaras);

		procOmegas.run();

		readGimp(rootOname, matsubaras, t);
	}

	void readGimp(PsimagLite::String filename,
	              const MatsubarasType& matsubaras,
	              DmrgType t)
	{
		std::ifstream fin(filename);
		if (!fin || !fin.good() || fin.bad())
			err("readGimp: Could not open " + filename + "\n");

		if (t == DmrgType::TYPE_0)
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

			if (t == DmrgType::TYPE_0)
				gimp_[ind] = ComplexType(val1, val2);
			else
				gimp_[ind] += ComplexType(val1, val2);

			++ind;

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
