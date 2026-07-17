#ifndef IMPURITYSOLVER_NEQ_TDMRG_H
#define IMPURITYSOLVER_NEQ_TDMRG_H

#include "CmdLineOptions.hh"
#include "DmrgRunner.h"
#include "ImpuritySolverNeqBase.h"
#include "KadanoffBaym.h"
#include "ParamsNeqDmftSolver.h"
#include <PsimagLite/PsimagLite.h>
#include <PsimagLite/Vector.h>
#include <cassert>
#include <cmath>
#include <complex>
#include <fstream>
#include <map>
#include <sstream>
#include <string>

namespace Dmft {

// Non-equilibrium impurity solver using DMRG++ time-dependent DMRG (tDMRG).
//
// Five-run approach using TargetingExpression with <P.last|P> gauge correction
// (from pr/dmrgpp-pvector-last-syntax).  This replaces the old TSPEvolveGroundState
// mechanism which did not fit the DMRG++ design.
//
// Run 1:  GS DMRG with H(U_i)  →  saves |GS_i⟩ to root_+gs.
// Run 2:  TargetingExpression restart, FiniteLoops flag=2 (keeps |GS_i⟩ as |gs⟩),
//         computes P0 = c'[0]*|gs⟩  →  saves to root_+particle.
// Run 3:  TargetingExpression tDMRG, restarts from root_+particle, flag=2:
//           P1 = TimeEvolve*|P0⟩ = e^{-iH_f t} c†|GS_i⟩
//           P2 = TimeEvolve*|gs⟩ = e^{-iH_f t} |GS_i⟩
//         In-situ: <P2|c|P1>   →  G^>(t,0)
//                  <P2.last|P2> →  consecutive GS overlap (gauge-phase correction)
// Run 4:  Same as Run 2 but with c[0]  →  saves to root_+hole.
// Run 5:  Same as Run 3 but restarts from root_+hole:
//         In-situ: <P1|c|P2>   →  G^<(t,0)
//                  <P2.last|P2> →  gauge-phase correction
//
// G^R(t,0) = G^>(t,0) - G^<(t,0)
//
// Validation: run with BandwidthFinal=BandwidthInitial and HubbardUFinal=HubbardU
// (no quench).  G^R(t,0) from tDMRG must agree with the ED solver to truncation error.
//
// Required input parameters:
//   TargetElectronsUp=, TargetElectronsDown=
//   RootOutputname=
//   InfiniteLoopKeptStates=
//   FiniteLoopsGs=        — loops for GS run (flag 0 OK)
//   FiniteLoopsTdmrg=     — loops for init/tDMRG runs (flag 2 enforced automatically;
//                           must contain enough sweeps for params_.nT time advances
//                           at TSPAdvanceEach sites each)
//   TSPTimeSteps=         — Krylov substep count (default 5, Gonzalo's preference)
//   TSPAdvanceEach=       — sites per time advance (use N-2 for N-site star; default 1)
template <typename ComplexOrRealType>
class ImpuritySolverNeqTdmrg : public ImpuritySolverNeqBase<ComplexOrRealType> {

public:

	using BaseType        = ImpuritySolverNeqBase<ComplexOrRealType>;
	using RealType        = typename BaseType::RealType;
	using ComplexType     = typename BaseType::ComplexType;
	using VectorRealType  = typename BaseType::VectorRealType;
	using KBType          = typename BaseType::KBType;
	using InputNgType     = typename BaseType::InputNgType;
	using ParamsNeqType   = ParamsNeqDmftSolver<ComplexOrRealType>;
	using DmrgRunnerType  = Dmrg::DmrgRunner<RealType>;
	using ApplicationType = PsimagLite::PsiApp;

	ImpuritySolverNeqTdmrg(const ParamsNeqType&            params,
	                       const ApplicationType&          app,
	                       typename InputNgType::Readable& io)
	    : params_(params)
	    , app_(app)
	    , io_(io)
	    , gimp_(params.nT,
	            params.eqParams.nMatsubaras,
	            params.dt,
	            params.eqParams.ficticiousBeta
	                / static_cast<RealType>(params.eqParams.nMatsubaras))
	{
		io.readline(nup_, "TargetElectronsUp=");
		io.readline(ndown_, "TargetElectronsDown=");
		io.readline(root_, "RootOutputname=");
		io.readline(infiniteLoops_, "InfiniteLoopKeptStates=");
		io.readline(finiteLoopsGs_, "FiniteLoopsGs=");

		try {
			io.readline(finiteLoopsTdmrg_, "FiniteLoopsTdmrg=");
		} catch (std::exception&) {
			finiteLoopsTdmrg_ = finiteLoopsGs_;
		}

		try {
			io.readline(tspTimeSteps_, "TSPTimeSteps=");
		} catch (std::exception&) {
			tspTimeSteps_ = 5;
		}

		try {
			io.readline(tspAdvanceEach_, "TSPAdvanceEach=");
		} catch (std::exception&) {
			tspAdvanceEach_ = 1;
		}
	}

	// Runs all five DMRG passes and fills the KB grid.
	void solve(const VectorRealType& bathParams) override
	{
		const SizeType nBath  = bathParams.size() / 2;
		const SizeType nsites = nBath + 1;

		VectorRealType hoppings(nBath), bathEps(nBath);
		for (SizeType i = 0; i < nBath; ++i) {
			hoppings[i] = bathParams[i];
			bathEps[i]  = bathParams[nBath + i];
		}

		VectorRealType potGS(nsites), potTdmrg(nsites);
		potGS[0]    = -RealType(0.5) * params_.uInitial;
		potTdmrg[0] = -RealType(0.5) * params_.uFinal;
		for (SizeType i = 0; i < nBath; ++i) {
			potGS[i + 1]    = bathEps[i];
			potTdmrg[i + 1] = bathEps[i];
		}

		// Run 1: GS DMRG for H(U_i)
		std::cout << "ImpuritySolverNeqTdmrg: GS DMRG U_i=" << params_.uInitial << "\n";
		{
			Dmrg::CmdLineOptions opts;
			opts.logfile = root_ + "gs.log";
			DmrgRunnerType runner(
			    app_,
			    buildGsInput(params_.uInitial, hoppings, potGS, nup_, ndown_, nsites),
			    opts);
			runner.doOneRun();
		}

		// Run 2: apply c'[0] to |GS_i⟩, save particle initial state
		std::cout << "ImpuritySolverNeqTdmrg: particle-init (c'|GS_i>)\n";
		{
			Dmrg::CmdLineOptions opts;
			opts.logfile = root_ + "particle_init.log";
			DmrgRunnerType runner(
			    app_,
			    buildParticleInitInput(
			        params_.uFinal, hoppings, potTdmrg, nup_, ndown_, nsites),
			    opts);
			runner.doOneRun();
		}

		// Run 3: particle tDMRG — G^>(t,0) and gauge
		const std::string tdmrgLog = root_ + "tdmrg.log";
		std::cout << "ImpuritySolverNeqTdmrg: particle tDMRG U_f=" << params_.uFinal
		          << " t_max=" << params_.tMax << "\n";
		{
			Dmrg::CmdLineOptions opts;
			opts.logfile              = tdmrgLog;
			opts.in_situ_measurements = "<P2|c|P1>,<P2.last|P2>";
			DmrgRunnerType runner(
			    app_,
			    buildTdmrgInput(
			        params_.uFinal, hoppings, potTdmrg, nup_, ndown_, nsites),
			    opts);
			runner.doOneRun();
		}

		// Run 4: apply c[0] to |GS_i⟩, save hole initial state
		std::cout << "ImpuritySolverNeqTdmrg: hole-init (c|GS_i>)\n";
		{
			Dmrg::CmdLineOptions opts;
			opts.logfile = root_ + "hole_init.log";
			DmrgRunnerType runner(
			    app_,
			    buildHoleInitInput(
			        params_.uFinal, hoppings, potTdmrg, nup_, ndown_, nsites),
			    opts);
			runner.doOneRun();
		}

		// Run 5: hole tDMRG — G^<(t,0) and gauge
		const std::string holeTdmrgLog = root_ + "tdmrg_hole.log";
		std::cout << "ImpuritySolverNeqTdmrg: hole tDMRG G^<(t,0)\n";
		{
			Dmrg::CmdLineOptions opts;
			opts.logfile              = holeTdmrgLog;
			opts.in_situ_measurements = "<P1|c|P2>,<P2.last|P2>";
			DmrgRunnerType runner(
			    app_,
			    buildHoleTdmrgInput(
			        params_.uFinal, hoppings, potTdmrg, nup_, ndown_, nsites),
			    opts);
			runner.doOneRun();
		}

		std::map<int, ComplexType> ggt0_at_step;
		parseTdmrgLog(tdmrgLog, ggt0_at_step);

		std::map<int, ComplexType> glt0_at_step;
		parseHoleTdmrgLog(holeTdmrgLog, glt0_at_step);

		fillKBGrid(ggt0_at_step, glt0_at_step);
	}

	void computeGimp(KBType& gimp, int n) const override
	{
		gimp.retarded(n, 0) = gimp_.retarded(n, 0);
		gimp.lesser(n, n)   = gimp_.lesser(n, n);
		gimp.lesser(n, 0)   = gimp_.lesser(n, 0);
	}

	const KBType& gimp() const override { return gimp_; }

private:

	// ---- Input construction ------------------------------------------------

	std::string buildGsInput(RealType              U,
	                         const VectorRealType& hoppings,
	                         const VectorRealType& potV,
	                         SizeType              nup,
	                         SizeType              ndown,
	                         SizeType              nsites) const
	{
		std::string s = "##Ainur1.0\n\n";
		s += geomHeader(nsites, U);
		s += "SolverOptions=twositedmrg,geometryallinsystem;\n";
		s += "Version=neqTdmrg;\n";
		s += "OutputFile=" + root_ + "gs;\n";
		s += "InfiniteLoopKeptStates=" + ttos(infiniteLoops_) + ";\n";
		s += "FiniteLoops=" + finiteLoopsGs_ + ";\n";
		s += "TargetElectronsUp=" + ttos(nup) + ";\n";
		s += "TargetElectronsDown=" + ttos(ndown) + ";\n";
		s += "dir0:Connectors=" + buildConnectorsStr(hoppings) + ";\n";
		s += "potentialV=" + buildPotentialVStr(potV) + ";\n";
		return s;
	}

	// Run 2: apply c'[0] to |GS_i⟩.
	// FiniteLoops flag=2 prevents DMRG from re-optimising |gs⟩ to |GS_f⟩,
	// keeping it as the WFT-transformed |GS_i⟩ throughout.
	std::string buildParticleInitInput(RealType              U_f,
	                                   const VectorRealType& hoppings,
	                                   const VectorRealType& potV,
	                                   SizeType              nup,
	                                   SizeType              ndown,
	                                   SizeType              nsites) const
	{
		std::string s = "##Ainur1.0\n\n";
		s += geomHeader(nsites, U_f);
		s += "SolverOptions=twositedmrg,geometryallinsystem,TargetingExpression,restart;\n";
		s += "Version=neqTdmrg;\n";
		s += "OutputFile=" + root_ + "particle;\n";
		s += "InfiniteLoopKeptStates=" + ttos(infiniteLoops_) + ";\n";
		s += "FiniteLoops=" + enforceFlag2(finiteLoopsGs_) + ";\n";
		s += "TargetElectronsUp=" + ttos(nup) + ";\n";
		s += "TargetElectronsDown=" + ttos(ndown) + ";\n";
		s += "dir0:Connectors=" + buildConnectorsStr(hoppings) + ";\n";
		s += "potentialV=" + buildPotentialVStr(potV) + ";\n";
		s += "RestartFilename=" + root_ + "gs;\n";
		s += "GsWeight=0.1;\n";
		s += "string P0=\"c'[0]*|gs>\";\n";
		return s;
	}

	// Run 3: particle tDMRG. P1=e^{-iH_f t}c†|GS_i⟩, P2=e^{-iH_f t}|GS_i⟩.
	// <P2|c|P1> = <GS_i(t)|c|Φ_+(t)⟩ → G^>(t,0) after multiplying by -i.
	// <P2.last|P2> = <GS_i(t-Δt)|GS_i(t)⟩ → N-sector gauge phase correction.
	std::string buildTdmrgInput(RealType              U_f,
	                            const VectorRealType& hoppings,
	                            const VectorRealType& potV,
	                            SizeType              nup,
	                            SizeType              ndown,
	                            SizeType              nsites) const
	{
		std::string s = "##Ainur1.0\n\n";
		s += geomHeader(nsites, U_f);
		s += "SolverOptions=twositedmrg,geometryallinsystem,TargetingExpression,restart,"
		     "usecomplex;\n";
		s += "Version=neqTdmrg;\n";
		s += "OutputFile=" + root_ + "tdmrg;\n";
		s += "InfiniteLoopKeptStates=" + ttos(infiniteLoops_) + ";\n";
		s += "FiniteLoops=" + enforceFlag2(finiteLoopsTdmrg_) + ";\n";
		s += "TargetElectronsUp=" + ttos(nup) + ";\n";
		s += "TargetElectronsDown=" + ttos(ndown) + ";\n";
		s += "dir0:Connectors=" + buildConnectorsStr(hoppings) + ";\n";
		s += "potentialV=" + buildPotentialVStr(potV) + ";\n";
		s += "RestartFilename=" + root_ + "particle;\n";
		s += "RestartMappingTvs=[0, -1, -1];\n";
		s += "GsWeight=0.1;\n";
		s += "string P0=|P0>;\n";
		s += "string P1=\"TimeEvolve{tau=" + ttos(params_.dt) + ",steps="
		    + ttos(tspTimeSteps_) + ",advanceEach=" + ttos(tspAdvanceEach_) + "}*|P0>\";\n";
		s += "string P2=\"TimeEvolve{tau=" + ttos(params_.dt) + ",steps="
		    + ttos(tspTimeSteps_) + ",advanceEach=" + ttos(tspAdvanceEach_) + "}*|gs>\";\n";
		return s;
	}

	// Run 4: apply c[0] to |GS_i⟩.  Mirror of Run 2 with annihilation operator.
	std::string buildHoleInitInput(RealType              U_f,
	                               const VectorRealType& hoppings,
	                               const VectorRealType& potV,
	                               SizeType              nup,
	                               SizeType              ndown,
	                               SizeType              nsites) const
	{
		std::string s = "##Ainur1.0\n\n";
		s += geomHeader(nsites, U_f);
		s += "SolverOptions=twositedmrg,geometryallinsystem,TargetingExpression,restart;\n";
		s += "Version=neqTdmrg;\n";
		s += "OutputFile=" + root_ + "hole;\n";
		s += "InfiniteLoopKeptStates=" + ttos(infiniteLoops_) + ";\n";
		s += "FiniteLoops=" + enforceFlag2(finiteLoopsGs_) + ";\n";
		s += "TargetElectronsUp=" + ttos(nup) + ";\n";
		s += "TargetElectronsDown=" + ttos(ndown) + ";\n";
		s += "dir0:Connectors=" + buildConnectorsStr(hoppings) + ";\n";
		s += "potentialV=" + buildPotentialVStr(potV) + ";\n";
		s += "RestartFilename=" + root_ + "gs;\n";
		s += "GsWeight=0.1;\n";
		s += "string P0=\"c[0]*|gs>\";\n";
		return s;
	}

	// Run 5: hole tDMRG.  P1=e^{-iH_f t}c|GS_i⟩, P2=e^{-iH_f t}|GS_i⟩.
	// <P1|c|P2> = <φ_h(t)|c|GS_i(t)⟩ → G^<(t,0) after multiplying by +i.
	std::string buildHoleTdmrgInput(RealType              U_f,
	                                const VectorRealType& hoppings,
	                                const VectorRealType& potV,
	                                SizeType              nup,
	                                SizeType              ndown,
	                                SizeType              nsites) const
	{
		std::string s = "##Ainur1.0\n\n";
		s += geomHeader(nsites, U_f);
		s += "SolverOptions=twositedmrg,geometryallinsystem,TargetingExpression,restart,"
		     "usecomplex;\n";
		s += "Version=neqTdmrg;\n";
		s += "OutputFile=" + root_ + "tdmrg_hole;\n";
		s += "InfiniteLoopKeptStates=" + ttos(infiniteLoops_) + ";\n";
		s += "FiniteLoops=" + enforceFlag2(finiteLoopsTdmrg_) + ";\n";
		s += "TargetElectronsUp=" + ttos(nup) + ";\n";
		s += "TargetElectronsDown=" + ttos(ndown) + ";\n";
		s += "dir0:Connectors=" + buildConnectorsStr(hoppings) + ";\n";
		s += "potentialV=" + buildPotentialVStr(potV) + ";\n";
		s += "RestartFilename=" + root_ + "hole;\n";
		s += "RestartMappingTvs=[0, -1, -1];\n";
		s += "GsWeight=0.1;\n";
		s += "string P0=|P0>;\n";
		s += "string P1=\"TimeEvolve{tau=" + ttos(params_.dt) + ",steps="
		    + ttos(tspTimeSteps_) + ",advanceEach=" + ttos(tspAdvanceEach_) + "}*|P0>\";\n";
		s += "string P2=\"TimeEvolve{tau=" + ttos(params_.dt) + ",steps="
		    + ttos(tspTimeSteps_) + ",advanceEach=" + ttos(tspAdvanceEach_) + "}*|gs>\";\n";
		return s;
	}

	// ---- Log parsing -------------------------------------------------------

	// Parse particle tDMRG log for <P2|c|P1> and <P2.last|P2>.
	// Applies gauge correction: ggt0 /= gauge_overlap.
	// Measurement format (site 0, from in_situ_measurements):
	//   "<site> (<re>,<im>) <time> <label>"
	void parseTdmrgLog(const std::string& logfile, std::map<int, ComplexType>& ggt0_at_step)
	{
		std::ifstream fin(logfile);
		if (!fin || !fin.good()) {
			err("ImpuritySolverNeqTdmrg: cannot open tDMRG log '" + logfile + "'\n");
			return;
		}

		std::map<int, ComplexType> gauge_at_step;

		std::string line;
		while (std::getline(fin, line)) {
			SizeType    site = 0;
			std::string valStr, label;
			RealType    t = 0;
			if (!parseMeasurementLine(line, site, valStr, t, label))
				continue;
			if (site != 0)
				continue;

			RealType re = 0, im = 0;
			if (!parseComplex(valStr, re, im))
				continue;

			const int n = static_cast<int>(std::round(t / params_.dt));
			if (n < 0 || n > static_cast<int>(params_.nT))
				continue;

			if (label == "<P2|c|P1>")
				ggt0_at_step[n] = ComplexType(re, im);
			else if (label == "<P2.last|P2>")
				gauge_at_step[n] = ComplexType(re, im);
		}

		if (ggt0_at_step.empty())
			std::cerr << "ImpuritySolverNeqTdmrg: WARNING: no G^> measurements in '"
			          << logfile << "'\n";

		applySignFlip(ggt0_at_step);

		// Divide by gauge: corrects N-sector phase drift between consecutive advances.
		for (auto& kv : ggt0_at_step) {
			auto it = gauge_at_step.find(kv.first);
			if (it != gauge_at_step.end() && std::abs(it->second) > RealType(1e-10))
				kv.second /= it->second;
		}

		applyGlobalPhase(ggt0_at_step);
	}

	// Parse hole tDMRG log for <P1|c|P2> and <P2.last|P2>.
	// Applies gauge correction: glt0 *= gauge_overlap.
	void parseHoleTdmrgLog(const std::string& logfile, std::map<int, ComplexType>& glt0_at_step)
	{
		std::ifstream fin(logfile);
		if (!fin || !fin.good()) {
			std::cerr << "ImpuritySolverNeqTdmrg: WARNING: cannot open hole log '"
			          << logfile << "'\n";
			return;
		}

		std::map<int, ComplexType> gauge_at_step;

		std::string line;
		while (std::getline(fin, line)) {
			SizeType    site = 0;
			std::string valStr, label;
			RealType    t = 0;
			if (!parseMeasurementLine(line, site, valStr, t, label))
				continue;
			if (site != 0)
				continue;

			RealType re = 0, im = 0;
			if (!parseComplex(valStr, re, im))
				continue;

			const int n = static_cast<int>(std::round(t / params_.dt));
			if (n < 0 || n > static_cast<int>(params_.nT))
				continue;

			if (label == "<P1|c|P2>")
				glt0_at_step[n] = ComplexType(re, im);
			else if (label == "<P2.last|P2>")
				gauge_at_step[n] = ComplexType(re, im);
		}

		if (glt0_at_step.empty()) {
			std::cerr << "ImpuritySolverNeqTdmrg: WARNING: no G^< measurements in '"
			          << logfile << "'\n";
			return;
		}

		applySignFlip(glt0_at_step);

		// Multiply by gauge: corrects N-sector phase for hole-run bra.
		for (auto& kv : glt0_at_step) {
			auto it = gauge_at_step.find(kv.first);
			if (it != gauge_at_step.end() && std::abs(it->second) > RealType(1e-10))
				kv.second *= it->second;
		}

		applyGlobalPhase(glt0_at_step);
	}

	void fillKBGrid(const std::map<int, ComplexType>& ggt0_at_step,
	                const std::map<int, ComplexType>& glt0_at_step)
	{
		const int         nT = static_cast<int>(params_.nT);
		const ComplexType pI = ComplexType(0, 1);
		const ComplexType mI = ComplexType(0, -1);

		for (int n = 0; n <= nT; ++n) {
			auto itG = ggt0_at_step.find(n);
			if (itG == ggt0_at_step.end())
				continue;
			const ComplexType ggt = mI * itG->second; // G^> = -i<P2|c|P1>

			auto itL = glt0_at_step.find(n);
			if (itL != glt0_at_step.end()) {
				const ComplexType glt = pI * itL->second; // G^< = +i<P1|c|P2>
				gimp_.lesser(n, 0)    = glt;
				gimp_.retarded(n, 0)  = ggt - glt;
			} else {
				gimp_.retarded(n, 0) = ggt;
			}
		}
	}

	// ---- Utilities ---------------------------------------------------------

	// Detect and correct sign flips from DMRG gauge jumps at sweep reversals.
	// A true sign flip satisfies |M(t)+M(t-dt)|² ≪ |M(t)-M(t-dt)|².
	static void applySignFlip(std::map<int, ComplexType>& m)
	{
		if (m.size() < 2)
			return;
		auto it   = m.begin();
		auto prev = it;
		++it;
		for (; it != m.end(); ++it, ++prev) {
			const ComplexType& a         = prev->second;
			const RealType     sum_norm  = std::norm(a + it->second);
			const RealType     diff_norm = std::norm(a - it->second);
			if (sum_norm < RealType(0.1) * diff_norm)
				it->second = -it->second;
		}
	}

	// Remove residual global phase: rotate so that the t=0 value is real.
	static void applyGlobalPhase(std::map<int, ComplexType>& m)
	{
		auto it0 = m.find(0);
		if (it0 == m.end())
			return;
		const ComplexType z0 = it0->second;
		if (std::abs(z0) < RealType(1e-10))
			return;
		const ComplexType phase_inv = ComplexType(std::abs(z0)) / z0;
		for (auto& kv : m)
			kv.second *= phase_inv;
	}

	// In-situ measurement log line: "<site> (<re>,<im>) <time> <label>"
	static bool parseMeasurementLine(const std::string& line,
	                                 SizeType&          site,
	                                 std::string&       valStr,
	                                 RealType&          t,
	                                 std::string&       label)
	{
		std::istringstream iss(line);
		if (!(iss >> site))
			return false;
		if (!(iss >> valStr))
			return false;
		if (!(iss >> t))
			return false;
		if (!(iss >> label))
			return false;
		return true;
	}

	static bool parseComplex(const std::string& s, RealType& re, RealType& im)
	{
		if (s.empty())
			return false;
		if (s[0] == '(') {
			const std::string body  = s.substr(1, s.size() > 2 ? s.size() - 2 : 0);
			auto              comma = body.find(',');
			if (comma == std::string::npos)
				return false;
			try {
				re = std::stod(body.substr(0, comma));
				im = std::stod(body.substr(comma + 1));
				return true;
			} catch (...) {
				return false;
			}
		}
		try {
			re = std::stod(s);
			im = RealType(0);
			return true;
		} catch (...) {
			return false;
		}
	}

	// Replace FiniteLoops flag 0 with flag 2 to prevent |gs⟩ re-optimisation.
	// The "2" (fast-WFT) flag keeps |gs⟩ = WFT-transformed |GS_i⟩, which is
	// required so that P2 = TimeEvolve*|gs⟩ evolves the pre-quench ground state.
	static std::string enforceFlag2(const std::string& loops)
	{
		std::string            result = loops;
		std::string::size_type pos    = 0;
		while ((pos = result.find(", 0]", pos)) != std::string::npos) {
			result.replace(pos, 4, ", 2]");
			pos += 4;
		}
		return result;
	}

	// Common geometry/model header block for all five runs.
	std::string geomHeader(SizeType nsites, RealType U) const
	{
		std::string s;
		s += "TotalNumberOfSites=" + ttos(nsites) + ";\n";
		s += "NumberOfTerms=1;\n";
		s += "DegreesOfFreedom=1;\n";
		s += "GeometryKind=star;\n";
		s += "GeometryOptions=none;\n";
		s += "hubbardU=" + buildHubbardUStr(U, nsites) + ";\n";
		s += "Model=HubbardOneBand;\n";
		return s;
	}

	static std::string buildHubbardUStr(RealType U, SizeType nsites)
	{
		std::string s = "[" + ttos(U);
		for (SizeType i = 1; i < nsites; ++i)
			s += ", 0.";
		return s + "]";
	}

	static std::string buildConnectorsStr(const VectorRealType& hoppings)
	{
		std::string s = "[";
		for (SizeType i = 0; i < hoppings.size(); ++i) {
			if (i > 0)
				s += ",";
			s += ttos(hoppings[i]);
		}
		return s + "]";
	}

	static std::string buildPotentialVStr(const VectorRealType& potV)
	{
		std::string inner;
		for (SizeType i = 0; i < potV.size(); ++i) {
			if (i > 0)
				inner += ",";
			inner += ttos(potV[i]);
		}
		return "[" + inner + "," + inner + "]";
	}

	// ---- Member data -------------------------------------------------------
	const ParamsNeqType&            params_;
	const ApplicationType&          app_;
	typename InputNgType::Readable& io_;
	KBType                          gimp_;
	SizeType                        nup_   = 0;
	SizeType                        ndown_ = 0;
	std::string                     root_;
	SizeType                        infiniteLoops_ = 0;
	std::string                     finiteLoopsGs_;
	std::string                     finiteLoopsTdmrg_;
	SizeType                        tspTimeSteps_   = 5;
	SizeType                        tspAdvanceEach_ = 1;
};

} // namespace Dmft
#endif // IMPURITYSOLVER_NEQ_TDMRG_H
