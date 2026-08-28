#include "CincuentaInputCheck.h"
#include "ParamsNeqDmftSolver.h"
#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_string.hpp>
#include <complex>
#include <string>

using ComplexType = std::complex<double>;
using ParamsType  = Dmft::ParamsNeqDmftSolver<ComplexType>;
using InputNgType = PsimagLite::InputNg<Dmft::CincuentaInputCheck>;
using Catch::Matchers::ContainsSubstring;

namespace {

const std::string kGbekInput = "##Ainur1.0\n\n"
                               "FicticiousBeta=16;\n"
                               "ChemicalPotential=1.;\n"
                               "Matsubaras=2;\n"
                               "LatticeGf=\"energy,semicircular,0\";\n"
                               "TargetElectronsUp=1;\n"
                               "TargetElectronsDown=0;\n"
                               "int ImpuritySite=0;\n"
                               "real HubbardU=2.;\n"
                               "HubbardUFinal=2.;\n"
                               "TmaxNeq=0.1;\n"
                               "NtNeq=1;\n"
                               "NeqDmftIter=1;\n"
                               "NeqBathRank=1;\n"
                               "NeqAtomicLimit=1;\n";

void requireValidGbekParams(const std::string& input)
{
	InputNgType::Writeable writable(Dmft::CincuentaInputCheck {}, input);
	InputNgType::Readable  readable(writable);
	Dmft::CincuentaInputCheck::rejectRemovedLabels(readable);
	ParamsType params(readable);
}

} // namespace

TEST_CASE("GBEK input requires a positive corrector count", "[GBEK][input]")
{
	auto input = kGbekInput;
	input.replace(
	    input.find("NeqDmftIter=1;"), std::string("NeqDmftIter=1;").size(), "NeqDmftIter=0;");
	CHECK_THROWS_WITH(requireValidGbekParams(input),
	                  ContainsSubstring("Non-equilibrium GBEK requires NeqDmftIter>0"));
}

TEST_CASE("GBEK input requires a positive bath rank", "[GBEK][input]")
{
	auto input = kGbekInput;
	input.replace(
	    input.find("NeqBathRank=1;"), std::string("NeqBathRank=1;").size(), "NeqBathRank=0;");
	CHECK_THROWS_WITH(requireValidGbekParams(input),
	                  ContainsSubstring("Non-equilibrium GBEK requires NeqBathRank>0"));
}

TEST_CASE("GBEK input rejects removed NeqSolver", "[GBEK][input]")
{
	const auto input = kGbekInput + "NeqSolver=\"exactdiag\";\n";
	CHECK_THROWS_WITH(requireValidGbekParams(input),
	                  ContainsSubstring("NeqSolver= was removed"));
}

TEST_CASE("GBEK input rejects removed NeqDmftTolerance", "[GBEK][input]")
{
	const auto input = kGbekInput + "NeqDmftTolerance=1e-8;\n";
	CHECK_THROWS_WITH(requireValidGbekParams(input),
	                  ContainsSubstring("NeqDmftTolerance= was removed"));
}
