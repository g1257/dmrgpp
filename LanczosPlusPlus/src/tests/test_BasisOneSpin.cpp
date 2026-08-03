#include "LanczosPlusPlus/src/Models/HubbardOneOrbital/BasisOneSpin.h"
#include <catch2/catch_test_macros.hpp>
#include <set>

using LanczosPlusPlus::BasisOneSpin;

// perfectIndex must be a bijection between a basis's stored words (index
// order) and [0, size): every stored word maps back to its own storage
// index, and no two stored words collide.
static void checkPerfectIndexIsBijection(const BasisOneSpin& basis)
{
	std::set<SizeType> seenIndices;
	for (SizeType i = 0; i < basis.size(); ++i) {
		const auto     word = basis[i];
		const SizeType idx  = basis.perfectIndex(word);
		INFO("i=" << i << " word=" << word << " perfectIndex=" << idx);
		CHECK(idx == i);
		CHECK(seenIndices.count(idx) == 0);
		seenIndices.insert(idx);
	}
}

TEST_CASE("BasisOneSpin perfectIndex is a bijection, small basis", "[BasisOneSpin]")
{
	BasisOneSpin basis(4, 2);
	checkPerfectIndexIsBijection(basis);
}

// Regression test for a real bug: BasisOneSpin's constructor calls
// LanczosGlobals::doCombinatorial(2*nsite+2), which resizes a
// process-lifetime static cache shared by every basis constructed anywhere
// in the same binary. A prior bug in that resize path left stale data
// behind whenever a basis with a smaller nsite was constructed after one
// with a larger nsite (see test_Combinatorial.cpp for the mechanism),
// silently corrupting perfectIndex for the smaller, later basis. This
// mirrors exactly how the bug surfaced in practice: as flaky Catch2 test
// failures that depended on which other tests (constructing differently
// sized systems) had already run in the same process.
TEST_CASE("BasisOneSpin perfectIndex stays correct after a larger basis was "
          "constructed first",
          "[BasisOneSpin][regression]")
{
	// Construct a deliberately larger system first, to perturb
	// LanczosGlobals' shared combinatorial cache to a larger size.
	{
		BasisOneSpin large(10, 5);
		checkPerfectIndexIsBijection(large);
	}

	// Now a smaller system: its combinatorial cache requirement (2*3+2=8)
	// is well below the larger one just used (2*10+2=22), so this
	// resizes the shared cache back down.
	BasisOneSpin small(3, 1);
	checkPerfectIndexIsBijection(small);

	// And smaller still, a couple more times, to catch any resize-order
	// sensitivity rather than just a single grow-then-shrink transition.
	BasisOneSpin smaller(2, 1);
	checkPerfectIndexIsBijection(smaller);

	BasisOneSpin tiny(1, 1);
	checkPerfectIndexIsBijection(tiny);
}
