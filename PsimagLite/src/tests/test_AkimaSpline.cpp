#include "AkimaSpline.h"
#include "AllocatorCpu.h"
#include <catch2/catch_approx.hpp>
#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers.hpp>
#include <catch2/matchers/catch_matchers_exception.hpp>
#include <catch2/matchers/catch_matchers_string.hpp>

namespace PsimagLite {
using Approx = Catch::Approx;

TEST_CASE("AkimaSpline_from_std_vector", "[psimaglite]")
{
	std::vector<double> x { 1.0, 2.0, 3.0, 4.0, 5.0 };
	std::vector<double> y { 1.0, 4.0, 9.0, 16.0, 25.0 };
	AkimaSpline akima_spline(x, y);
	CHECK(akima_spline(1.0) == Approx(1.0));
	CHECK(akima_spline(2.0) == Approx(4.0));
	CHECK(akima_spline(3.0) == Approx(9.0));
	CHECK(akima_spline(3.5) == Approx(12.25));
	CHECK(akima_spline(4.0) == Approx(16.0));
	CHECK(akima_spline(5.0) == Approx(25.0));
	CHECK_THROWS_MATCHES(akima_spline(25.0),
	                     RuntimeError,
	                     MessageMatches(Catch::Matchers::StartsWith("Akima: out of range")));
}
}
