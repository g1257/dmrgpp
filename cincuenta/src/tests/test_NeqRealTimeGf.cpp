#include "NeqRealTimeGf.h"
#include <PsimagLite/PsimagLite.h>
#include <catch2/catch_test_macros.hpp>
#include <cstdio>
#include <fstream>
#include <string>

using Dmft::NeqRealTimeGf;

TEST_CASE("NeqRealTimeGf owns only the real-time GBEK components", "[NeqRealTimeGf]")
{
	NeqRealTimeGf<double> gf(4, 0.25);
	CHECK(gf.nT() == 4);
	CHECK(gf.dt() == 0.25);
	CHECK(SizeType(gf.retarded.rows()) == 5);
	CHECK(SizeType(gf.retarded.cols()) == 5);
	CHECK(SizeType(gf.lesser.rows()) == 5);
	CHECK(SizeType(gf.lesser.cols()) == 5);
}

TEST_CASE("NeqRealTimeGf dump writes the retained real-time format", "[NeqRealTimeGf]")
{
	NeqRealTimeGf<double> gf(1, 0.5);
	gf.retarded(0, 0) = { 0.0, -1.0 };
	gf.retarded(1, 0) = { 0.25, -0.75 };
	gf.retarded(1, 1) = { 0.0, -1.0 };
	gf.lesser(0, 0)  = { 0.0, 0.5 };
	gf.lesser(0, 1)  = { 0.1, 0.2 };
	gf.lesser(1, 0)  = { -0.1, 0.2 };
	gf.lesser(1, 1)  = { 0.0, 0.5 };

	const std::string prefix = "test-neq-real-time-gf";
	gf.dump(prefix);

	std::ifstream retarded(prefix + "-retarded");
	std::ifstream lesser(prefix + "-lesser");
	std::string line;
	REQUIRE(std::getline(retarded, line));
	CHECK(line == "0.0000000000 0.0000000000 0.0000000000 -1.0000000000");
	REQUIRE(std::getline(lesser, line));
	CHECK(line == "0.0000000000 0.0000000000 0.0000000000 0.5000000000");

	std::remove((prefix + "-retarded").c_str());
	std::remove((prefix + "-lesser").c_str());
}
