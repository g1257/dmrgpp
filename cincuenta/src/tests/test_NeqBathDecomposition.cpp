#include "KadanoffBaym.h"
#include "NeqBathDecomposition.h"
#include <catch2/catch_approx.hpp>
#include <catch2/catch_test_macros.hpp>
#include <cmath>
#include <complex>
#include <vector>

using Dmft::KadanoffBaym;
using Dmft::NeqBathDecomposition;

// ── Helpers ───────────────────────────────────────────────────────────────

namespace {

// NeqBathDecomposition with a single dummy bath site (V=0, ε=0).
// This makes deltaMinusLesser = 0, so -iΔ⁺<(n,j) = -i * delta.lesser(n,j).
static NeqBathDecomposition<double> makeDecomp(SizeType rank, SizeType nT)
{
	PsimagLite::Vector<double>::Type bathParams = { 0.0, 0.0 };
	return NeqBathDecomposition<double>(rank,
	                                    /*beta=*/10.0,
	                                    /*mu=*/0.0,
	                                    bathParams,
	                                    nT,
	                                    /*nTau=*/4,
	                                    /*dt=*/0.2,
	                                    /*dtau=*/0.5);
}

// Fill delta.lesser so that -iΔ⁺<(n,j) = target[n][j].
// With zero first bath: iDeltaPlusLesser = i * delta.lesser,
// so setting delta.lesser(n,j) = i * target gives -iDeltaPlusLesser = target. ✓
static KadanoffBaym<double> makeDelta(SizeType nT, const std::vector<std::vector<double>>& target)
{
	const std::complex<double> I(0, 1);
	KadanoffBaym<double>       delta(nT, /*nTau=*/4, /*dt=*/0.2, /*dtau=*/0.5);
	for (SizeType n = 0; n <= nT; ++n)
		for (SizeType j = 0; j <= nT; ++j)
			delta.lesser(n, j) = I * target[n][j];
	return delta;
}

// Call update for n = 0..nT in sequence.
static void
runUpdates(NeqBathDecomposition<double>& decomp, const KadanoffBaym<double>& delta, SizeType nT)
{
	for (int n = 0; n <= static_cast<int>(nT); ++n)
		decomp.update(n, delta);
}

// Max |Σ_p V(n,p) conj(V(j,p)) − target[n][j]| over all n,j.
static double reconstructionError(NeqBathDecomposition<double>&           decomp,
                                  SizeType                                nT,
                                  const std::vector<std::vector<double>>& target)
{
	double maxErr = 0.0;
	for (SizeType n = 0; n <= nT; ++n) {
		for (SizeType j = 0; j <= nT; ++j) {
			std::complex<double> val(0);
			for (SizeType p = 0; p < decomp.rank(); ++p)
				val += decomp.Vplus(static_cast<int>(n), static_cast<int>(p))
				    * std::conj(
				           decomp.Vplus(static_cast<int>(j), static_cast<int>(p)));
			maxErr = std::max(maxErr, std::abs(val.real() - target[n][j]));
		}
	}
	return maxErr;
}

// Construct a physically-motivated PSD target matrix of rank r:
//   target[n][j] = Σ_{k=0}^{r-1} v[k][n] * v[k][j]
// with v[k][0] = 0 (degenerate at n=0, matching DMFT physics: Δ⁺<(0,j) = 0).
static std::vector<std::vector<double>> makePsdTarget(SizeType nT, int rank)
{
	const SizeType                   N = nT + 1;
	std::vector<std::vector<double>> v(rank, std::vector<double>(N, 0.0));
	for (int r = 0; r < rank; ++r) {
		const double decay = 0.3 + 0.2 * r;
		const double freq  = 0.5 * (r + 1);
		for (SizeType n = 1; n < N; ++n)
			v[r][n] = std::exp(-decay * n) * std::sin(freq * n);
	}

	std::vector<std::vector<double>> mat(N, std::vector<double>(N, 0.0));
	for (SizeType n = 0; n < N; ++n)
		for (SizeType j = 0; j < N; ++j)
			for (int r = 0; r < rank; ++r)
				mat[n][j] += v[r][n] * v[r][j];
	return mat;
}

} // namespace

// ── Seeding / structural tests ─────────────────────────────────────────────

TEST_CASE("n=0 is always skipped — V_(0,p) stays zero", "[NeqBathDecomposition][seeding]")
{
	const SizeType nT = 6;
	for (SizeType L : { 1u, 2u, 3u }) {
		auto target = makePsdTarget(nT, static_cast<int>(L));
		auto delta  = makeDelta(nT, target);
		auto decomp = makeDecomp(L, nT);
		runUpdates(decomp, delta, nT);

		for (SizeType p = 0; p < L; ++p)
			CHECK(std::abs(decomp.Vplus(0, static_cast<int>(p))) < 1e-15);
	}
}

TEST_CASE("L=2: column 0 is alive after fix", "[NeqBathDecomposition][seeding]")
{
	// Column 0 is seeded at n=1 (pivot row 1, not 0). With the old code
	// V_(0,0)≈0 caused a guard-to-zero and column 0 stayed dead (L_eff=1).
	const SizeType nT     = 6;
	auto           target = makePsdTarget(nT, 2);
	auto           delta  = makeDelta(nT, target);
	auto           decomp = makeDecomp(2, nT);
	runUpdates(decomp, delta, nT);

	CHECK(std::abs(decomp.Vplus(1, 0)) > 1e-6); // V_(1,0) non-zero
	CHECK(std::abs(decomp.Vplus(1, 1)) < 1e-15); // col 1 not yet seeded at n=1
	CHECK(std::abs(decomp.Vplus(2, 1)) > 1e-6); // V_(2,1) seeded at n=2
}

TEST_CASE("L=3: all three columns alive after seeding phase", "[NeqBathDecomposition][seeding]")
{
	const SizeType nT     = 8;
	auto           target = makePsdTarget(nT, 3);
	auto           delta  = makeDelta(nT, target);
	auto           decomp = makeDecomp(3, nT);
	runUpdates(decomp, delta, nT);

	CHECK(std::abs(decomp.Vplus(1, 0)) > 1e-6); // col 0 seeded at n=1
	CHECK(std::abs(decomp.Vplus(2, 1)) > 1e-6); // col 1 seeded at n=2
	CHECK(std::abs(decomp.Vplus(3, 2)) > 1e-6); // col 2 seeded at n=3

	// Columns not yet seeded at their seed step must be zero
	CHECK(std::abs(decomp.Vplus(1, 1)) < 1e-15);
	CHECK(std::abs(decomp.Vplus(1, 2)) < 1e-15);
	CHECK(std::abs(decomp.Vplus(2, 2)) < 1e-15);
}

// ── Reconstruction accuracy tests ─────────────────────────────────────────

TEST_CASE("L=1 reconstructs a rank-1 matrix exactly", "[NeqBathDecomposition][reconstruction]")
{
	const SizeType nT     = 10;
	auto           target = makePsdTarget(nT, 1);
	auto           delta  = makeDelta(nT, target);
	auto           decomp = makeDecomp(1, nT);
	runUpdates(decomp, delta, nT);

	// A rank-1 Cholesky of a rank-1 matrix should be exact up to floating point.
	CHECK(reconstructionError(decomp, nT, target) < 1e-10);
}

TEST_CASE("L=2 reconstructs a rank-2 matrix accurately", "[NeqBathDecomposition][reconstruction]")
{
	const SizeType nT     = 10;
	auto           target = makePsdTarget(nT, 2);
	auto           delta  = makeDelta(nT, target);
	auto           decomp = makeDecomp(2, nT);
	runUpdates(decomp, delta, nT);

	// The optimal update minimizes the residual in the column space; for an
	// exactly rank-2 input the reconstruction should be essentially exact.
	CHECK(reconstructionError(decomp, nT, target) < 1e-8);
}

TEST_CASE("L=1 has larger error than L=2 on a rank-2 matrix",
          "[NeqBathDecomposition][reconstruction]")
{
	// This is the regression test that would have caught the L_eff=1 bug:
	// if rank-2 runs secretly deliver rank-1 accuracy, this test fails.
	const SizeType nT     = 10;
	auto           target = makePsdTarget(nT, 2);
	auto           delta  = makeDelta(nT, target);

	auto decomp1 = makeDecomp(1, nT);
	runUpdates(decomp1, delta, nT);
	const double err1 = reconstructionError(decomp1, nT, target);

	auto decomp2 = makeDecomp(2, nT);
	runUpdates(decomp2, delta, nT);
	const double err2 = reconstructionError(decomp2, nT, target);

	// L=2 must be strictly better than L=1.
	CHECK(err2 < err1);
	// And L=1's residual must be measurable (not numerically zero).
	CHECK(err1 > 1e-4);
}

TEST_CASE("L=3 reconstructs a rank-3 matrix accurately", "[NeqBathDecomposition][reconstruction]")
{
	const SizeType nT     = 12;
	auto           target = makePsdTarget(nT, 3);
	auto           delta  = makeDelta(nT, target);
	auto           decomp = makeDecomp(3, nT);
	runUpdates(decomp, delta, nT);

	CHECK(reconstructionError(decomp, nT, target) < 1e-8);
}

// ── Incisive tests: full-rank target, many steps past the seeding phase ────
//
// All tests above use makePsdTarget(), which is EXACTLY rank-r and L == r
// throughout. That construction cannot distinguish a correct "optimal
// update" step (GBEK Eq. 62-63) from either of two bugs this test was built
// to catch -- once L == rank, the residual is zero regardless of how the
// update step is implemented, so both bugs below went undetected by the
// pre-existing tests:
//
//  1. Fixed reference window: the least-squares reference set was
//     hardcoded to the first L seeding rows forever, instead of growing to
//     ALL previously-determined rows every step. Causes near-total
//     collapse of the reconstructed diagonal within a handful of steps
//     past L.
//  2. Missing diagonal constraint: after fixing (1), the update step
//     solved only the linear part of Eq. 63 (Q^H Q q = Q^H a), ignoring
//     the q^H q ~ d term that the paper's actual objective
//     F(q) = 2||Qq-a||^2 + |q^H q - d|^2 requires jointly. Undershoots the
//     diagonal by a smaller amount, growing with n.
//
// The target below is an Ornstein-Uhlenbeck / exponential kernel,
// 0.5*exp(-decay*|n-j|*dt) (PSD by Bochner's theorem, and genuinely
// full-rank -- not reducible to any finite rank exactly), run for 30 steps
// at rank L=3, i.e. 27 steps into the "optimal update" phase, far more than
// any existing test here exercises.
//
// Expected values were generated independently in Python, directly from
// the paper's own Eq. 56-63 (NOT ported from this C++ file), via:
//
//   cd cincuenta/TestSuite/gbek_reference
//   uv run --with numpy python3 -c "
//   import numpy as np
//   from gbek_cholesky import cholesky_causal
//   dt, decay, N = 0.2, 0.15, 31
//   target = np.zeros((N, N))
//   for n in range(N):
//       for j in range(N):
//           target[n, j] = 0.5 * np.exp(-decay * abs(n - j) * dt)
//   target[0, :] = 0.0
//   target[:, 0] = 0.0
//   V = cholesky_causal(target.astype(complex), 3)
//   diag = (V @ V.conj().T).real.diagonal()
//   for n in [1, 5, 10, 15, 20, 25, 30]:
//       print(n, repr(diag[n]))
//   "
//
// The same file's cholesky_causal_buggy_fixed_window() (bug 1) and
// cholesky_causal_linear_only() (bug 2 only, bug 1 already fixed) reproduce
// the two intermediate stages of the bug for regression/contrast -- see
// that file's module docstring for the full before/after story, and its
// own __main__ block for a side-by-side table of target vs. fixed vs.
// linear-only vs. buggy.
namespace {

static std::vector<std::vector<double>> makeFullRankTarget(SizeType nT, double dt, double decay)
{
	const SizeType                   N = nT + 1;
	std::vector<std::vector<double>> mat(N, std::vector<double>(N, 0.0));
	for (SizeType n = 0; n < N; ++n)
		for (SizeType j = 0; j < N; ++j)
			mat[n][j] = 0.5
			    * std::exp(-decay * std::abs(static_cast<int>(n) - static_cast<int>(j))
			               * dt);
	for (SizeType j = 0; j < N; ++j)
		mat[0][j] = 0.0;
	for (SizeType n = 0; n < N; ++n)
		mat[n][0] = 0.0;
	return mat;
}

} // namespace

TEST_CASE("Full-rank target: L=3 reconstruction matches the Python reference "
          "(GBEK Eq. 62-63), not the fixed-window bug",
          "[NeqBathDecomposition][reconstruction][incisive]")
{
	const SizeType nT     = 30;
	const double   dt     = 0.2;
	const double   decay  = 0.15;
	auto           target = makeFullRankTarget(nT, dt, decay);
	auto           delta  = makeDelta(nT, target);
	auto           decomp = makeDecomp(3, nT);
	runUpdates(decomp, delta, nT);

	auto diagAt = [&](SizeType n)
	{
		std::complex<double> val(0);
		for (SizeType p = 0; p < decomp.rank(); ++p)
			val += decomp.Vplus(static_cast<int>(n), static_cast<int>(p))
			    * std::conj(decomp.Vplus(static_cast<int>(n), static_cast<int>(p)));
		return val.real();
	};

	// Values generated by gbek_cholesky.py's cholesky_causal(), per the
	// command in the comment above -- this is the correct-algorithm answer
	// (both bugs fixed), not what either the fixed-window bug (1) or the
	// linear-only bug (2) produced (see that file's __main__ table).
	CHECK(diagAt(1) == Catch::Approx(0.5000000000000001).epsilon(1e-6));
	CHECK(diagAt(5) == Catch::Approx(0.49021304302869706).epsilon(1e-6));
	CHECK(diagAt(10) == Catch::Approx(0.4771665196067349).epsilon(1e-6));
	CHECK(diagAt(15) == Catch::Approx(0.4637717267240376).epsilon(1e-6));
	CHECK(diagAt(20) == Catch::Approx(0.45078669197353066).epsilon(1e-6));
	CHECK(diagAt(25) == Catch::Approx(0.437773730548557).epsilon(1e-6));
	CHECK(diagAt(30) == Catch::Approx(0.42485714581863027).epsilon(1e-6));

	// Discriminating bound: the buggy fixed-window algorithm (bug 1) gives
	// 0.099 at n=30 (min diag ~0.099 over n=1..30) on this exact target --
	// this assertion alone would have failed under that implementation,
	// independent of the exact-value checks above. It would NOT, by
	// itself, have caught bug 2 (linear-only gives 0.419 at n=30, still
	// above 0.35) -- that one needs the tight-tolerance exact-value checks.
	CHECK(diagAt(30) > 0.35);
}

// ── Complex-phase target: catches the swapped-conjugate bug ────────────────
//
// Every target above (makePsdTarget, makeFullRankTarget) is REAL-valued, so
// V ends up real too (a real symmetric PSD matrix has a real Cholesky
// factor), making conj(V) == V identically. That made all of them blind to
// a genuine third bug: choleskyOptimalUpdate's Gram matrix and RHS
// conjugated the WRONG factor (QtQ(p,pp) = V(k,p)*conj(V(k,pp)) instead of
// conj(V(k,p))*V(k,pp); Qta with no conjugate on V at all). This is
// mathematically well-defined (QtQ stays Hermitian either way) so it never
// crashed or produced NaN -- it silently solved a different linear system
// whenever V's columns carry genuine complex phase, exactly the situation
// in any real GBEK run (the atomic-limit and near-atomic targets both
// oscillate in time). Found 2026-07-08 by dumping cincuenta's actual
// online-computed V_(n,p) row-by-row and comparing against
// gbek_cholesky.py::cholesky_causal() fed the identical (Q,a,d) inputs --
// see cincuenta/TestSuite/gbek_reference/compare_V_rows.py.
namespace {

static std::vector<std::vector<std::complex<double>>>
makeComplexPhaseTarget(SizeType nT, double dt, double decay, double omega)
{
	const SizeType                                 N = nT + 1;
	std::vector<std::vector<std::complex<double>>> mat(
	    N, std::vector<std::complex<double>>(N, std::complex<double>(0)));
	const std::complex<double> I(0, 1);
	for (SizeType n = 1; n < N; ++n) {
		for (SizeType j = 1; j < N; ++j) {
			const double tn = n * dt;
			const double tj = j * dt;
			mat[n][j]       = 0.5 * std::exp(-decay * std::abs(tn - tj))
			    * std::exp(I * omega * (tn - tj));
		}
	}
	return mat;
}

static KadanoffBaym<double>
makeComplexDelta(SizeType nT, const std::vector<std::vector<std::complex<double>>>& target)
{
	const std::complex<double> I(0, 1);
	KadanoffBaym<double>       delta(nT, /*nTau=*/4, /*dt=*/0.2, /*dtau=*/0.5);
	for (SizeType n = 0; n <= nT; ++n)
		for (SizeType j = 0; j <= nT; ++j)
			delta.lesser(n, j) = I * target[n][j];
	return delta;
}

} // namespace

TEST_CASE("Complex-phase target: L=3 reconstruction matches the Python reference "
          "(catches the swapped-conjugate bug in QtQ/Qta)",
          "[NeqBathDecomposition][reconstruction][incisive]")
{
	const SizeType nT     = 30;
	const double   dt     = 0.2;
	const double   decay  = 0.4;
	const double   omega  = 1.2;
	auto           target = makeComplexPhaseTarget(nT, dt, decay, omega);
	auto           delta  = makeComplexDelta(nT, target);
	auto           decomp = makeDecomp(3, nT);
	runUpdates(decomp, delta, nT);

	// Values generated by gbek_cholesky.py's cholesky_causal() (correct
	// conjugation convention) on the identical target -- command:
	//   python3 -c "from gbek_cholesky import cholesky_causal; ..." with
	//   the target built as above (decay=0.4, omega=1.2, dt=0.2).
	//
	// Regenerated 2026-07-10 (see choleskyOptimalUpdate's bug-3 comment):
	// the previous set of hardcoded values here (from 2026-07-09) was
	// itself generated from a gbek_cholesky.py that had the SAME swapped-
	// conjugate bug this test's name warns about -- both that Python
	// snapshot and this file's pre-2026-07-10 C++ used the same wrong
	// convention, so they agreed with each other while both being wrong.
	// These values are from the actually-correct convention, verified
	// three independent ways (exact rank-1 atomic-limit reconstruction,
	// lambda-scaling convergence to an analytic O(v^2) benchmark, and
	// landing within a few percent of GBEK Fig. 8's quoted peak). The
	// OLD (2026-07-09) hardcoded values differ dramatically from these at
	// every row past the seeding phase (n>3), e.g. at n=30 the two
	// conventions differ by max|diff|~0.2 -- nowhere close to this
	// tolerance, exactly the signature bug 3 produces.
	auto checkRow = [&](SizeType             n,
	                    std::complex<double> v0,
	                    std::complex<double> v1,
	                    std::complex<double> v2)
	{
		CHECK(decomp.Vplus(static_cast<int>(n), 0).real()
		      == Catch::Approx(v0.real()).epsilon(1e-6));
		CHECK(decomp.Vplus(static_cast<int>(n), 0).imag()
		      == Catch::Approx(v0.imag()).epsilon(1e-6));
		CHECK(decomp.Vplus(static_cast<int>(n), 1).real()
		      == Catch::Approx(v1.real()).epsilon(1e-6));
		CHECK(decomp.Vplus(static_cast<int>(n), 1).imag()
		      == Catch::Approx(v1.imag()).epsilon(1e-6));
		CHECK(decomp.Vplus(static_cast<int>(n), 2).real()
		      == Catch::Approx(v2.real()).epsilon(1e-6));
		CHECK(decomp.Vplus(static_cast<int>(n), 2).imag()
		      == Catch::Approx(v2.imag()).epsilon(1e-6));
	};

	checkRow(4,
	         { 0.4208631578224682, 0.36912556596604357 },
	         { 0.16775833438485793, 0.08733680808421918 },
	         { 0.35551055925393726, 0.08699937184079912 });
	checkRow(10,
	         { -0.25222285556232754, 0.3773515074833713 },
	         { 0.12378637598989402, -0.3399544834009459 },
	         { -0.034962857490067334, 0.3188882144219492 });
	checkRow(20,
	         { -0.06106508225086797, -0.39761179617963466 },
	         { 0.15398515251391467, 0.3720793521156658 },
	         { 0.13205440664592072, 0.1802101423453549 });
	checkRow(30,
	         { 0.2776917711887221, 0.22309940380070978 },
	         { -0.08896792402662543, -0.04153868424300278 },
	         { -0.41272431535331866, -0.08229556453955512 });
}
