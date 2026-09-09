#include "CincuentaInputCheck.h"
#include "DmftSolver.h"
#include <catch2/catch_approx.hpp>
#include <catch2/catch_test_macros.hpp>
#include <cmath>
#include <cstdio>
#include <fstream>
#include <string>

using ComplexType = std::complex<double>;
using SolverType  = Dmft::DmftSolver<ComplexType>;
using ParamsType  = SolverType::ParamsDmftSolverType;
using InputNgType = PsimagLite::InputNg<Dmft::CincuentaInputCheck>;

namespace {

const std::string config = "##Ainur1.0\n\n"
                           "FicticiousBeta=10.;\n"
                           "ChemicalPotential=1.;\n"
                           "Matsubaras=3;\n"
                           "LatticeGf=\"energy,semicircular,4\";\n"
                           "NumberOfBathPoints=1;\n"
                           "DmftNumberOfIterations=1;\n"
                           "DmftTolerance=1e-3;\n"
                           "ImpuritySolver=\"exactdiag\";\n"
                           "FitOptions=particleholesymmetric;\n"
                           "MinParamsDelta=0.01;\n"
                           "MinParamsDelta2=0.01;\n"
                           "MinParamsTolerance=1e-4;\n"
                           "MinParamsMaxIter=100;\n"
                           "MinParamsVerbose=0;\n"
                           "TargetElectronsUp=1;\n"
                           "TargetElectronsDown=1;\n"
                           "int ImpuritySite=0;\n"
                           "real HubbardU=2.;\n"
                           "RootOutputname=\"testEquilibriumInitialData\";\n"
                           "InfiniteLoopKeptStates=20;\n"
                           "matrix FiniteLoopsGs=[[@auto, 20, 0],[@auto, 20, 0]];\n"
                           "real OmegaBegin=-4.;\n"
                           "integer OmegaTotal=7;\n"
                           "real OmegaStep=0.2;\n"
                           "real OmegaDelta=0.2;\n"
                           "integer TridiagSteps=100;\n"
                           "real TridiagEps=1e-6;\n"
                           "TruncationTolerance=\"1e-6,20\";\n"
                           "CorrectionVectorEta=0.;\n"
                           "GsWeight=0.1;\n"
                           "matrix FiniteLoopsOmega=[[@auto, 20, 2],[@auto, 20, 2]];\n";

} // namespace

TEST_CASE("DmftSolver snapshots equilibrium Matsubara data before real solve",
          "[DmftSolver][EquilibriumInitialData]")
{
	InputNgType::Writeable           ioWriteable(Dmft::CincuentaInputCheck {}, config);
	InputNgType::Readable            io(ioWriteable);
	ParamsType                       params(io);
	SolverType::FitType::InitResults initResults(io);

	int                argc    = 1;
	char               name[]  = "test_DmftSolverEquilibriumInitialData";
	char*              argv[]  = { name, nullptr };
	char**             argvPtr = argv;
	PsimagLite::PsiApp app("test_DmftSolverEquilibriumInitialData", &argc, &argvPtr, 1);

	SolverType solver(params, initResults, app, io);
	solver.selfConsistencyLoop();
	const auto& handoff = solver.equilibriumInitialData();

	// OmegaTotal=7 makes the final real-frequency gimp cardinality different
	// from the 2*Matsubaras=6 equilibrium data retained by the handoff.
	REQUIRE(handoff.gimpMatsubara.size() == 2 * params.nMatsubaras);
	REQUIRE(handoff.matsubaraFrequencies.size() == handoff.gimpMatsubara.size());
	REQUIRE(handoff.bathParameters.size() == solver.bathResult().size());

	for (SizeType i = 0; i < handoff.bathParameters.size(); ++i)
		CHECK(handoff.bathParameters[i] == Catch::Approx(solver.bathResult()[i]));

	for (SizeType i = 0; i < handoff.matsubaraFrequencies.size(); ++i) {
		const double      expected = (i < params.nMatsubaras)
		         ? -M_PI * (2 * (params.nMatsubaras - i) - 1) / params.ficticiousBeta
		         : M_PI * (2 * (i - params.nMatsubaras) + 1) / params.ficticiousBeta;
		const ComplexType value    = handoff.gimpMatsubara[i];
		const ComplexType mate
		    = handoff.gimpMatsubara[handoff.gimpMatsubara.size() - 1 - i];
		CHECK(handoff.matsubaraFrequencies[i] == Catch::Approx(expected));
		CHECK(std::isfinite(std::real(value)));
		CHECK(std::isfinite(std::imag(value)));
		CHECK(std::abs(value) > 1e-12);
		CHECK(std::real(value) == Catch::Approx(std::real(mate)));
		CHECK(std::imag(value) == Catch::Approx(-std::imag(mate)));
	}

	const std::string filename = "equilibrium-gimp-matsubara";
	handoff.writeMatsubara(filename);
	std::ifstream input(filename);
	REQUIRE(input);
	SizeType rows      = 0;
	double   frequency = 0;
	double   real      = 0;
	double   imag      = 0;
	while (input >> frequency >> real >> imag) {
		REQUIRE(rows < handoff.matsubaraFrequencies.size());
		CHECK(frequency == Catch::Approx(handoff.matsubaraFrequencies[rows]));
		CHECK(real == Catch::Approx(std::real(handoff.gimpMatsubara[rows])));
		CHECK(imag == Catch::Approx(std::imag(handoff.gimpMatsubara[rows])));
		++rows;
	}
	CHECK(input.eof());
	CHECK(rows == handoff.matsubaraFrequencies.size());
	std::remove(filename.c_str());
}
