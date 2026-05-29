#ifndef IMPURITYSOLVER_DMRG_H
#define IMPURITYSOLVER_DMRG_H

#include "CmdLineOptions.hh"
#include "DmrgRunner.h"
#include "Geometry/Star.h"
#include "ImpuritySolverBase.h"
#include "InputCheck.h"
#include "LanczosSolver.h"
#include "ManyOmegas.h"
#include "Matsubaras.h"
#include "ModelParams.h"
#include "ParamsDmftSolver.h"
#include "ProcOmegas.h"
#include "PsiBase64.h"
#include "PsimagLite.h"
#include "Vector.h"

namespace Dmft {

template <typename ComplexOrRealType>
class ImpuritySolverDmrg : public ImpuritySolverBase<ComplexOrRealType> {

	enum class DmrgType
	{
		TYPE_0,
		TYPE_1
	};

public:

	using BaseType               = ImpuritySolverBase<ComplexOrRealType>;
	using ParamsDmftSolverType   = ParamsDmftSolver<ComplexOrRealType>;
	using RealType               = typename PsimagLite::Real<ComplexOrRealType>::Type;
	using ComplexType            = std::complex<RealType>;
	using VectorRealType         = typename PsimagLite::Vector<RealType>::Type;
	using VectorComplexType      = typename PsimagLite::Vector<ComplexType>::Type;
	using DmrgRunnerType         = Dmrg::DmrgRunner<RealType>;
	using ApplicationType        = PsimagLite::PsiApp;
	using MatsubarasType         = PsimagLite::Matsubaras<RealType>;
	using RealFrequencyRangeType = PsimagLite::RealFrequencyRange<RealType>;
	using ModelParamsType        = typename BaseType::ModelParamsType;
	using InputNgType            = typename BaseType::InputNgType;

	ImpuritySolverDmrg(const ParamsDmftSolverType&     params,
	                   const ApplicationType&          app,
	                   typename InputNgType::Readable& io)
	    : BaseType(params.ficticiousBeta, params.nMatsubaras, io)
	    , params_(params)
	    , app_(app)
	    , io_(io)
	    , freq_enum_(PsimagLite::FreqEnum::MATSUBARA)
	{ }

	// bathParams[0-nBath-1] ==> V ==> hoppings impurity --> bath
	// bathParams[nBath-...] ==> energies on each bath site
	void solve(const VectorRealType& bathParams, PsimagLite::FreqEnum freq_enum, SizeType iter)
	{
		ModelParamsType model_params(bathParams, io_);
		SizeType        mpiRank = PsimagLite::MPI::commRank(PsimagLite::MPI::COMM_WORLD);

		if (mpiRank == 0) {
			PsimagLite::String data2 = BaseType::createGsInput(model_params, io_);
			// PsimagLite::String insitu = "<gs|nup|gs>";

			Dmrg::CmdLineOptions cmdline_options;
			cmdline_options.in_situ_measurements = "<gs|nup|gs>";
			cmdline_options.logfile              = "-";

			DmrgRunnerType runner(app_, data2, cmdline_options);
			runner.doOneRun();
		}

		PsimagLite::MPI::barrier(PsimagLite::MPI::COMM_WORLD);

		PsimagLite::String data3 = createOmegaInput(model_params, freq_enum);

		SizeType impurity_site = model_params.impuritySite();

		doType(DmrgType::TYPE_0, data3, impurity_site, freq_enum);

		doType(DmrgType::TYPE_1, data3, impurity_site, freq_enum);

		BaseType::scale(gimp_);

		freq_enum_ = freq_enum;

		PsimagLite::MPI::barrier(PsimagLite::MPI::COMM_WORLD);
	}

	const VectorComplexType& gimp() const { return gimp_; }

	PsimagLite::FreqEnum freqEnum() const { return freq_enum_; }

private:

	std::string createOmegaInput(const ModelParamsType& model_params,
	                             PsimagLite::FreqEnum   freq_enum) const
	{
		std::string s = BaseType::commonInputString(
		    model_params, io_, BaseType::GsOrOmegaEnum::OMEGA);

		std::string root;
		io_.readline(root, "RootOutputname=");
		s += "RestartFilename=" + root + "gs;\n";

		try {
			SizeType tsteps = 0;
			io_.readline(tsteps, "TridiagSteps=");
			s += "int TridiagSteps=" + ttos(tsteps) + ";\n";
		} catch (std::exception&) { }

		try {
			RealType teps = 0;
			io_.readline(teps, "TridiagEps=");
			s += "real TridiagEps=" + ttos(teps) + ";\n";
		} catch (std::exception&) { }

		try {
			std::string tt;
			io_.readline(tt, "TruncationTolerance=");
			s += "TruncationTolerance=" + tt + ";\n";
		} catch (std::exception&) { }

		s += "CorrectionA=0;\n";

		if (freq_enum == PsimagLite::FreqEnum::MATSUBARA) {
			s += "CorrectionVectorFreqType=Matsubara;\n";
		} else {
			s += "CorrectionVectorFreqType=Real;\n";
		}

		RealType eta = 0;
		if (freq_enum == PsimagLite::FreqEnum::REAL) {
			io_.readline(eta, "OmegaDelta=");
		} else {
			io_.readline(eta, "CorrectionVectorEta=");
		}

		s += "CorrectionVectorEta=" + ttos(eta) + ";\n";

		s += "CorrectionVectorAlgorithm=Krylov;\n";
		s += "Orbitals=1;\nGsWeight=0.1;\n";

		try {
			RealType gsw = 0;
			io_.readline(gsw, "GsWeight=");
			s += "GsWeight=" + ttos(gsw) + ";\n";
		} catch (std::exception&) { }

		std::string data = BaseType::addBathParams(s, model_params);

		std::ofstream tout("testout.ain");
		tout << data;
		tout.close();
		return data;
	}

	/*
	 * This business of communicating by strings is far from ideal,
	 * but DMRG++ doesn't have an internal API right now
	 */
	static std::string
	addTypeAndObs(DmrgType t, SizeType impurity_site, PsimagLite::String data)
	{
		const std::string obsTc = (t == DmrgType::TYPE_0) ? "c'" : "c";
		return data + addType(t) + addObs(impurity_site, obsTc);
	}

	/* The type added is because there are two terms that need to be computed
	 * by separate runs and the terms differ by a sign
	 */
	static std::string addType(DmrgType t)
	{

		const SizeType tt = (t == DmrgType::TYPE_0) ? 0 : 1;
		return "DynamicDmrgType=" + ttos(tt) + ";\n";
	}

	/*
	 * The observable is just one for the Green function: c (the destruction operator)
	 * and we apply at the impurity. Now, if the impurity is site 0 (a border site)
	 * we need to add a trigger at site 1, that is, and operator at site 1: the identity,
	 * that is we multiply by the identity in this case.
	 */
	static std::string addObs(SizeType impurity_site, const std::string& obs)
	{
		std::string str = "TSPProductOrSum=product;\n";
		SizeType    ind = 0;
		if (impurity_site > 0) {
			str += "TSPSites=[" + ttos(impurity_site)
			    + "];\n"
			      " TSPLoops=[0]\n";
			ind = 0;
		} else {
			str += "TSPSites=[1, 0];\n"
			       "TSPLoops=[0, 0];\n"
			       "string TSPOp0:TSPOperator=expression;\n"
			       "string TSPOp0:OperatorExpression=identity;\n"
			       "string TSPOp1:TSPOperator=expression;\n";
			ind = 1;
		}

		str += "string TSPOp" + ttos(ind) + ":OperatorExpression=" + obs + ";\n";
		return str;
	}

	void doType(DmrgType             t,
	            PsimagLite::String   data,
	            SizeType             impurity_site,
	            PsimagLite::FreqEnum freq_enum)
	{
		std::string obs   = (t == DmrgType::TYPE_0) ? "c" : "c'";
		std::string data2 = addTypeAndObs(t, impurity_site, data);

		runOmegas(data2, obs, freq_enum);

		procOmegas(data2, t, freq_enum);
	}

	void procOmegas(const std::string& data2, DmrgType t, PsimagLite::FreqEnum freq_enum)
	{
		SizeType mpiRank = PsimagLite::MPI::commRank(PsimagLite::MPI::COMM_WORLD);

		if (mpiRank != 0)
			return;

		const std::string rootIname   = "input";
		const std::string rootOname   = "OUTPUT";
		const bool        skipFourier = true;

		Dmrg::InputCheck                                          inputCheck;
		typename PsimagLite::InputNg<Dmrg::InputCheck>::Writeable ioW(inputCheck, data2);
		typename PsimagLite::InputNg<Dmrg::InputCheck>::Readable  io(ioW);

		SizeType total = 0;
		if (freq_enum == PsimagLite::FreqEnum::MATSUBARA) {
			total = this->matsubaras().total();
			Dmrg::ProcOmegas<RealType, MatsubarasType> procOmegas(io,
			                                                      params_.precision,
			                                                      skipFourier,
			                                                      rootIname,
			                                                      rootOname,
			                                                      this->matsubaras());

			procOmegas.run();
		} else {
			total = this->realFreqRange().total();
			Dmrg::ProcOmegas<RealType, RealFrequencyRangeType> procOmegas(
			    io,
			    params_.precision,
			    skipFourier,
			    rootIname,
			    rootOname,
			    this->realFreqRange());

			procOmegas.run();
		}

		readGimp(rootOname, total, t);
	}

	void runOmegas(const std::string&   data2,
	               const std::string&   obs,
	               PsimagLite::FreqEnum freq_enum) const
	{
		const bool               dryrun   = false;
		const PsimagLite::String rootname = "dmftDynamics";
		Dmrg::CmdLineOptions     cmdline_options;
		cmdline_options.in_situ_measurements = "<gs|" + obs + "|P2>,<gs|" + obs + "|P3>";

		if (freq_enum == PsimagLite::FreqEnum::MATSUBARA) {
			Dmrg::ManyOmegas<RealType, MatsubarasType> manyOmegas(
			    data2, this->matsubaras(), app_);
			manyOmegas.run(dryrun, rootname, cmdline_options);
		} else {
			Dmrg::ManyOmegas<RealType, RealFrequencyRangeType> manyOmegas(
			    data2, this->realFreqRange(), app_);
			manyOmegas.run(dryrun, rootname, cmdline_options);
		}
	}

	void readGimp(PsimagLite::String filename, SizeType total, DmrgType t)
	{
		std::ifstream fin(filename);
		if (!fin || !fin.good() || fin.bad())
			err("readGimp: Could not open " + filename + "\n");

		gimp_.resize(total);

		SizeType ind = 0;
		while (!fin.eof()) {
			RealType val = 0;
			fin >> val;
			SizeType n = 0;
			fin >> n;

			bool     centerSeen = false;
			RealType val1       = 0;
			RealType val2       = 0;

			for (SizeType i = 0; i <= params_.center_site; ++i) {
				SizeType site = 0;
				fin >> site;

				fin >> val2;

				fin >> val1;

				if (site == params_.center_site) {
					centerSeen = true;
					break;
				}
			}

			if (!centerSeen)
				err("Internal error: center " + ttos(params_.center_site)
				    + " not seen, freq id = " + ttos(ind) + "\n");

			assert(ind < gimp_.size());
			if (t == DmrgType::TYPE_0) {
				gimp_[ind] = ComplexType(-val1, -val2);
			} else {
				gimp_[ind] += ComplexType(-val1, -val2);
			}

			++ind;

			if (ind >= gimp_.size())
				break;

			if (n == 1)
				continue;

			const SizeType tmp = params_.center_site + 1;
			n -= tmp;

			for (SizeType i = 0; i < 3 * n; ++i)
				fin >> val1;
		}

		if (ind < gimp_.size())
			err("readGimp: Not all values computed\n");
	}

	const ParamsDmftSolverType&     params_;
	const ApplicationType&          app_;
	typename InputNgType::Readable& io_;
	VectorComplexType               gimp_;
	PsimagLite::FreqEnum            freq_enum_;
};
}
#endif // IMPURITYSOLVER_DMRG_H
