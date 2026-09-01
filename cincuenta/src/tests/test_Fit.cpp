#include "Fit.h"
#include "LatticeGf.h"
#include <algorithm>
#include <catch2/catch_test_macros.hpp>
#include <cmath>
#include <complex>
#include <stdexcept>

namespace {

/**
 * @brief Minimal readable configuration adapter for constructing `MinParams` in unit tests.
 *
 * `MinParams` accepts any duck-typed input object that provides `readline(value, label)`.
 * This fixture supplies the fitting parameters used by `inputCincuenta0.ain` without
 * requiring a complete Ainur input file or an `InputNg` parser. Each overload handles
 * the value type requested by `MinParams`; unsupported labels throw because `MinParams`
 * interprets that exception as "parameter not provided" and retains its default. In
 * particular, rejecting `FitMethod=` selects the default
 * conjugate-gradient method.
 * Note there is both a FitMethod and a FitOption in play for fit and
 * only FitMethod is part of DMFT::Fit's expectation for a minParams
 * passed to it.
 * Future tests of `Dmft::Fit` can reuse this pattern, adding supported labels only when
 * the behavior under test requires values that differ from the production defaults.
 */
struct FitConfig {
	void readline(double& value, const PsimagLite::String& label)
	{
		if (label == "MinParamsDelta=")
			value = 0.1;
		else if (label == "MinParamsDelta2=")
			value = 0.01;
		else if (label == "MinParamsTolerance=")
			value = 1e-3;
		else
			throw std::runtime_error("unknown floating-point fit parameter");
	}

	void readline(SizeType& value, const PsimagLite::String& label)
	{
		if (label == "MinParamsMaxIter=")
			value = 10000;
		else
			throw std::runtime_error("unknown integer fit parameter");
	}

	void readline(int& value, const PsimagLite::String& label)
	{
		if (label == "MinParamsVerbose=")
			value = 0;
		else
			throw std::runtime_error("unknown integer fit parameter");
	}

	void readline(PsimagLite::String&, const PsimagLite::String&)
	{
		throw std::runtime_error("use the default fit method");
	}
};

} // namespace

TEST_CASE("centered PH fit produces a finite noncollapsed bath", "[Fit][particle-hole]")
{
	using ComplexType             = std::complex<double>;
	using FunctionOfFrequencyType = Dmft::FunctionOfFrequency<ComplexType>;
	using FitType                 = Dmft::Fit<ComplexType>;

	constexpr SizeType           nBath = 5;
	FunctionOfFrequencyType      sigma(/*fictitiousBeta=*/50.0, /*nMatsubara=*/400);
	Dmft::LatticeGf<ComplexType> latticeGf(sigma, /*mu=*/0.0, "energy,semicircular,4");
	latticeGf.update();

	FitConfig                  config;
	FitType::MinParamsType     minParams(config);
	const FitType::InitResults initialBath(/*ra=*/1.0, /*rb=*/0.0, {}, /*reset=*/false);
	FitType                    fit(nBath, minParams, initialBath);
	fit.fit(latticeGf.g0(), /*mu=*/0.0, FitType::Options::PARTICLE_HOLE_SYMM);

	const auto& bath = fit.result();
	REQUIRE(bath.size() == 2 * nBath);
	for (const double value : bath)
		REQUIRE(std::isfinite(value));

	const auto hopping = [&bath](SizeType i) { return bath[i]; };
	const auto energy  = [&bath](SizeType i) { return bath[nBath + i]; };

	CHECK(hopping(3) == hopping(0));
	CHECK(hopping(4) == hopping(1));
	CHECK(energy(2) == 0.0);
	CHECK(energy(3) == -energy(0));
	CHECK(energy(4) == -energy(1));

	double maxHybridization = 0.0;
	for (SizeType i = 0; i < nBath; ++i)
		maxHybridization = std::max(maxHybridization, std::abs(hopping(i)));
	CHECK(maxHybridization > 1e-6);
}
