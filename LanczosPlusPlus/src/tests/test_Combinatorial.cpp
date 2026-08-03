#include "LanczosPlusPlus/src/Engine/Combinatorial.hh"
#include <catch2/catch_test_macros.hpp>

using LanczosPlusPlus::Combinatorial;

// Independent reference, not sharing Combinatorial's own algorithm.
static unsigned long expectedBinomial(unsigned n, unsigned k)
{
	if (k > n)
		return 0;
	unsigned long result = 1;
	for (unsigned i = 0; i < k; ++i)
		result = result * (n - i) / (i + 1);
	return result;
}

TEST_CASE("Combinatorial matches C(n,k) for a fresh table", "[Combinatorial]")
{
	Combinatorial c(9);
	for (SizeType n = 0; n < 9; ++n)
		for (SizeType k = 0; k < 9; ++k)
			CHECK(c(n, k) == expectedBinomial(n, k));
}

// Regression test for a real bug: Combinatorial::resize used to leave stale
// data behind when shrinking to a size it hadn't been at before, because it
// called Matrix::resize(m,m,0) (not element-preserving across a row-count
// change) instead of assigning a fresh Matrix. C(n,k) for k>n is never
// explicitly written by doCombinatorial's fill loop (only k<=n is
// meaningful), so those cells silently depended on whatever was left behind
// by a previous, differently-shaped table -- e.g. shrinking from 12x12 to
// 8x8 left the new (0,1) holding the old (8,0)'s value (1) instead of the
// correct 0. This broke LanczosPlusPlus::BasisHubbardLanczos::perfectIndex
// for any Fock word with its lowest set bit at position 0 (i.e. routinely),
// whenever two differently-sized bases were constructed in the same process
// -- exactly what happens across Catch2 test cases sharing this one binary.
TEST_CASE("Combinatorial stays correct across repeated resizes to different sizes",
          "[Combinatorial][regression]")
{
	Combinatorial c;
	c.resize(1);
	c.resize(12); // grow well past any size used elsewhere in this test
	c.resize(8); // shrink to a size that was NOT the original size

	for (SizeType n = 0; n < 8; ++n)
		for (SizeType k = 0; k < 8; ++k)
			CHECK(c(n, k) == expectedBinomial(n, k));

	// The specific cell the real bug corrupted.
	CHECK(c(0, 1) == 0);

	// Grow again, shrink to yet another size never seen before.
	c.resize(20);
	c.resize(5);
	for (SizeType n = 0; n < 5; ++n)
		for (SizeType k = 0; k < 5; ++k)
			CHECK(c(n, k) == expectedBinomial(n, k));
}
