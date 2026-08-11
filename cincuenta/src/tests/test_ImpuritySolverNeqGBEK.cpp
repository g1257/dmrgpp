#include "CincuentaInputCheck.h"
#include "ImpuritySolverNeqGBEK.h"
#include "TestMatrixUtils.h"
#include <catch2/catch_test_macros.hpp>
#include <cmath>

// ── Types ────────────────────────────────────────────────────────────────────

using ComplexType   = std::complex<double>;
using RealType      = double;
using SolverType    = Dmft::ImpuritySolverNeqGBEK<ComplexType>;
using ParamsType    = Dmft::ParamsNeqDmftSolver<ComplexType>;
using InputNgType   = PsimagLite::InputNg<Dmft::CincuentaInputCheck>;
using VectorComplex = typename PsimagLite::Vector<ComplexType>::Type;
using CrsMatrixType = PsimagLite::CrsMatrix<ComplexType>;
using WordType      = LanczosPlusPlus::LanczosGlobals::WordType;
using Acc           = Dmft::GBEKTestAccessor;

// ── Ainur config ─────────────────────────────────────────────────────────────
// Minimal GBEK config: 1 first-bath site, 2 time steps, rank-1 second bath.
// Equilibrium solver is ExactDiag; GBEK path invoked for neq propagation.
// DMRG fields are included because ParamsDmftSolver requires them even when
// the neq solver is GBEK.
static const std::string kConfig = "##Ainur1.0\n\n"
                                   "FicticiousBeta=10;\n"
                                   "ChemicalPotential=0.;\n"
                                   "Matsubaras=20;\n"
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
                                   "HubbardUFinal=2.;\n"
                                   "RootOutputname=\"testGBEK\";\n"
                                   "InfiniteLoopKeptStates=20;\n"
                                   "matrix FiniteLoopsGs=[[@auto, 20, 0],[@auto, 20, 0]];\n"
                                   "real OmegaBegin=-4.;\n"
                                   "integer OmegaTotal=40;\n"
                                   "real OmegaStep=0.2;\n"
                                   "real OmegaDelta=0.2;\n"
                                   "integer TridiagSteps=100;\n"
                                   "real TridiagEps=1e-6;\n"
                                   "TruncationTolerance=\"1e-6,20\";\n"
                                   "CorrectionVectorEta=0.;\n"
                                   "GsWeight=0.1;\n"
                                   "matrix FiniteLoopsOmega=[[@auto, 20, 2],[@auto, 20, 2]];\n"
                                   "TmaxNeq=0.2;\n"
                                   "NtNeq=2;\n"
                                   "NeqDmftIter=1;\n"
                                   "NeqDmftTolerance=1e-4;\n"
                                   "NeqSolver=\"gbek\";\n"
                                   "NeqBathRank=1;\n"
                                   "BandwidthFinal=3.;\n";

// ── Shared fixture ────────────────────────────────────────────────────────────
// Constructed once (static) so the full solve() runs only on first access.
// Members are declared in construction order so lifetimes are nested correctly:
//   ioW -> io -> params -> solver.
struct GBEKFixture {
	InputNgType::Writeable ioW;
	InputNgType::Readable  io;
	ParamsType             params;
	SolverType             solver;

	GBEKFixture()
	    : ioW(Dmft::CincuentaInputCheck {}, kConfig)
	    , io(ioW)
	    , params(io)
	    , solver(params, io)
	{
		// V=0.5, epsilon=0 for the single first-bath site.
		solver.solve({ 0.5, 0.0 });
	}
};

static GBEKFixture& fixture()
{
	static GBEKFixture f;
	return f;
}

// ── Helpers ───────────────────────────────────────────────────────────────────

// Returns a test vMid with bathRank entries using known complex values.
static std::vector<ComplexType> makeTestVMid(SizeType bathRank)
{
	std::vector<ComplexType> v(bathRank);
	for (SizeType p = 0; p < bathRank; ++p)
		v[p] = ComplexType(0.3 + 0.1 * p, 0.2 - 0.05 * p);
	return v;
}

// Compare CSR SpMV against applyHext for a single (v, vMid) combination.
// Returns max|hv_csr - hv_ref|.
static RealType spMVDiff(const SolverType&               solver,
                         const VectorComplex&            v,
                         const std::vector<ComplexType>& vMid,
                         const CrsMatrixType&            csr,
                         const std::vector<WordType>&    upW,
                         const std::vector<WordType>&    dnW,
                         SizeType                        dim1)
{
	const SizeType dim = upW.size() * dnW.size();
	VectorComplex  hv_csr(dim, ComplexType(0));
	VectorComplex  hv_ref(dim, ComplexType(0));

	csr.matrixVectorProduct(hv_csr, v);
	Acc::applyHext(solver, v, hv_ref, vMid, upW, dnW, dim1);

	RealType maxDiff = 0;
	for (SizeType i = 0; i < dim; ++i)
		maxDiff = std::max(maxDiff, std::abs(hv_csr[i] - hv_ref[i]));
	return maxDiff;
}

// ── TC-1: applyHext is Hermitian ─────────────────────────────────────────────
// Validates the reference implementation before using it to check anything else.
// If this fails, the physics was broken before the CSR existed.
TEST_CASE("applyHext is Hermitian for N-1 and N+1 sectors", "[GBEK][applyHext]")
{
	auto&          f = fixture();
	auto&          s = f.solver;
	const SizeType r = Acc::bathRank(s);

	SECTION("N-1 sector, zero vMid (fixed part only)")
	{
		std::vector<ComplexType> zeroVMid(r, ComplexType(0));
		auto                     H = buildDenseFromApplyHext(
                    s, zeroVMid, Acc::upWordsNm1(s), Acc::dnWordsNm1(s), Acc::dim1Nm1(s));
		assertHermitian(H);
	}

	SECTION("N-1 sector, non-zero vMid")
	{
		auto vMid = makeTestVMid(r);
		auto H    = buildDenseFromApplyHext(
                    s, vMid, Acc::upWordsNm1(s), Acc::dnWordsNm1(s), Acc::dim1Nm1(s));
		assertHermitian(H);
	}

	SECTION("N+1 sector, non-zero vMid")
	{
		auto vMid = makeTestVMid(r);
		auto H    = buildDenseFromApplyHext(
                    s, vMid, Acc::upWordsNp1(s), Acc::dnWordsNp1(s), Acc::dim1Np1(s));
		assertHermitian(H);
	}
}

// ── TC-2: buildHextCSR produces a Hermitian matrix ───────────────────────────
// If TC-1 passes but TC-2 fails the bug is in buildHextCSR, not applyHext.
TEST_CASE("buildHextCSR produces Hermitian matrix for both sectors", "[GBEK][buildHextCSR]")
{
	auto& f = fixture();
	auto& s = f.solver;

	SECTION("N-1 sector after updateCSR with test vMid")
	{
		auto vMid = makeTestVMid(Acc::bathRank(s));
		Acc::updateCSR(s, Acc::csrNm1Mut(s), Acc::varNm1(s), vMid);
		auto H = buildDenseFromCSR(Acc::csrNm1(s));
		assertHermitian(H);
	}

	SECTION("N+1 sector after updateCSR with test vMid")
	{
		auto vMid = makeTestVMid(Acc::bathRank(s));
		Acc::updateCSR(s, Acc::csrNp1Mut(s), Acc::varNp1(s), vMid);
		auto H = buildDenseFromCSR(Acc::csrNp1(s));
		assertHermitian(H);
	}
}

// ── TC-3: updateCSR writes variable slots and preserves fixed slots ───────────
TEST_CASE("updateCSR modifies only variable slots", "[GBEK][updateCSR]")
{
	auto&          f    = fixture();
	auto&          s    = f.solver;
	const SizeType r    = Acc::bathRank(s);
	auto           vMid = makeTestVMid(r);

	// Work on N-1 sector.
	auto&          csr    = Acc::csrNm1Mut(s);
	const auto&    varVec = Acc::varNm1(s);
	const SizeType nnz    = csr.nonZeros();

	// Snapshot all values before update.
	std::vector<ComplexType> before(nnz);
	for (SizeType i = 0; i < nnz; ++i)
		before[i] = csr.getValue(i);

	Acc::updateCSR(s, csr, varVec, vMid);

	// Build set of variable indices for O(1) membership test.
	std::vector<bool> isVar(nnz, false);
	for (const auto& ve : varVec)
		isVar[ve.nnzIdx] = true;

	for (SizeType i = 0; i < nnz; ++i) {
		if (isVar[i]) {
			// Expected value computed from VarEntry.
			// (We check it via the aggregate SpMV test; here we just verify it
			// changed.) A variable slot set from a non-zero vMid should differ from its
			// zero build-time value unless the VarEntry happens to produce zero. Accept
			// either "changed" (normal) or "correctly zero" (edge case).
			(void)i; // checked in TC-4
		} else {
			// Fixed slots must be unchanged.
			INFO("Fixed slot " << i << " was modified by updateCSR");
			CHECK(csr.getValue(i) == before[i]);
		}
	}

	// Verify each variable slot has the value that the VarEntry formula predicts.
	for (const auto& ve : varVec) {
		const ComplexType Vp       = ve.isConj ? std::conj(vMid[ve.p]) : vMid[ve.p];
		const ComplexType expected = Vp * ComplexType(static_cast<RealType>(ve.sign));
		INFO("VarEntry nnzIdx=" << ve.nnzIdx << " p=" << ve.p << " isConj=" << ve.isConj
		                        << " sign=" << ve.sign);
		CHECK(std::abs(csr.getValue(ve.nnzIdx) - expected) < 1e-15);
	}
}

// ── TC-4: SpMV matches applyHext for multiple (v, vMid) combinations ─────────
TEST_CASE("CSR SpMV agrees with applyHext for multiple vectors and vMid", "[GBEK][SpMV][applyHext]")
{
	auto&          f = fixture();
	auto&          s = f.solver;
	const SizeType r = Acc::bathRank(s);

	// Work with the N-1 sector.
	const auto&    upW    = Acc::upWordsNm1(s);
	const auto&    dnW    = Acc::dnWordsNm1(s);
	const SizeType d      = upW.size() * dnW.size();
	const SizeType d1     = Acc::dim1Nm1(s);
	auto&          csr    = Acc::csrNm1Mut(s);
	const auto&    varVec = Acc::varNm1(s);

	REQUIRE(d > 0);

	// Three vMid cases.
	std::vector<std::vector<ComplexType>> vMidCases = {
	    std::vector<ComplexType>(r, ComplexType(0)),     // zero
	    makeTestVMid(r),                                  // complex test values
	    [r]() {                                           // purely imaginary
	        std::vector<ComplexType> v(r);
	        for (SizeType p = 0; p < r; ++p)
	            v[p] = ComplexType(0, 0.7 + 0.1 * p);
	        return v;
	    }()
	};

	// Three test vectors.
	VectorComplex e0(d, ComplexType(0));
	e0[0] = ComplexType(1);
	VectorComplex emid(d, ComplexType(0));
	emid[d / 2] = ComplexType(1);
	// Simple normalised linear combination.
	VectorComplex esum(d, ComplexType(0));
	if (d >= 2) {
		esum[0]     = ComplexType(1) / std::sqrt(2.0);
		esum[d - 1] = ComplexType(1) / std::sqrt(2.0);
	} else {
		esum[0] = ComplexType(1);
	}

	const std::vector<const VectorComplex*> vecs = { &e0, &emid, &esum };

	for (SizeType vi = 0; vi < vMidCases.size(); ++vi) {
		Acc::updateCSR(s, csr, varVec, vMidCases[vi]);
		for (SizeType ti = 0; ti < vecs.size(); ++ti) {
			const RealType diff
			    = spMVDiff(s, *vecs[ti], vMidCases[vi], csr, upW, dnW, d1);
			INFO("vMid case " << vi << ", test vec " << ti
			                  << ": max|CSR*v - applyHext(v)| = " << diff);
			CHECK(diff < 1e-12);
		}
	}
}

// ── TC-5: krylovExpmvCSR agrees with krylovExpmv ────────────────────────────
TEST_CASE("krylovExpmvCSR agrees with krylovExpmv to 1e-10", "[GBEK][Krylov]")
{
	auto&          f    = fixture();
	auto&          s    = f.solver;
	const SizeType r    = Acc::bathRank(s);
	auto           vMid = makeTestVMid(r);

	const auto&    upW    = Acc::upWordsNm1(s);
	const auto&    dnW    = Acc::dnWordsNm1(s);
	const SizeType d      = upW.size() * dnW.size();
	const SizeType d1     = Acc::dim1Nm1(s);
	auto&          csr    = Acc::csrNm1Mut(s);
	const auto&    varVec = Acc::varNm1(s);

	// dt from TmaxNeq/NtNeq = 0.2/2 = 0.1
	const RealType dt = 0.1;

	Acc::updateCSR(s, csr, varVec, vMid);

	auto checkPair = [&](const VectorComplex& psi, const char* label)
	{
		auto result_ref = Acc::krylovExpmv(s, psi, vMid, upW, dnW, d1, dt);
		auto result_csr = Acc::krylovExpmvCSR(s, psi, Acc::csrNm1(s), dt);

		REQUIRE(result_ref.size() == d);
		REQUIRE(result_csr.size() == d);

		RealType maxDiff = 0;
		for (SizeType i = 0; i < d; ++i)
			maxDiff = std::max(maxDiff, std::abs(result_csr[i] - result_ref[i]));

		INFO(label << ": max|krylovCSR - krylovRef| = " << maxDiff);
		CHECK(maxDiff < 1e-10);
	};

	// Test vector 1: standard basis e_0.
	VectorComplex e0(d, ComplexType(0));
	e0[0] = ComplexType(1);
	checkPair(e0, "e_0");

	// Test vector 2: bStates[0] (the seeded c_{imp,up}|GS> state).
	const auto& hist = Acc::bStates(s);
	if (!hist.empty() && hist[0].size() == d) {
		checkPair(hist[0], "bStates[0]");
	}
}

// ── TC-6: krylovExpmvCSR preserves norm ──────────────────────────────────────
// A non-Hermitian H causes norm drift; this catches wrong-sign JW entries
// or numerical errors in the Krylov recurrence.
TEST_CASE("krylovExpmvCSR preserves vector norm to 1e-10", "[GBEK][Krylov][norm]")
{
	auto&          f    = fixture();
	auto&          s    = f.solver;
	const SizeType r    = Acc::bathRank(s);
	auto           vMid = makeTestVMid(r);

	const auto&    upW    = Acc::upWordsNm1(s);
	const auto&    dnW    = Acc::dnWordsNm1(s);
	const SizeType d      = upW.size() * dnW.size();
	auto&          csr    = Acc::csrNm1Mut(s);
	const auto&    varVec = Acc::varNm1(s);

	const RealType dt = 0.1;

	Acc::updateCSR(s, csr, varVec, vMid);

	// Unit vector e_0.
	VectorComplex e0(d, ComplexType(0));
	e0[0] = ComplexType(1);

	auto           result  = Acc::krylovExpmvCSR(s, e0, Acc::csrNm1(s), dt);
	const RealType normOut = vecNorm(result);
	INFO("||krylovExpmvCSR(e_0)|| = " << normOut);
	CHECK(std::abs(normOut - 1.0) < 1e-10);

	// Also check N+1 sector to exercise both paths.
	const auto&    upWp = Acc::upWordsNp1(s);
	const auto&    dnWp = Acc::dnWordsNp1(s);
	const SizeType dp   = upWp.size() * dnWp.size();
	auto&          csrP = Acc::csrNp1Mut(s);
	const auto&    varP = Acc::varNp1(s);

	Acc::updateCSR(s, csrP, varP, vMid);

	VectorComplex ep(dp, ComplexType(0));
	ep[0] = ComplexType(1);

	auto           resultP = Acc::krylovExpmvCSR(s, ep, Acc::csrNp1(s), dt);
	const RealType normP   = vecNorm(resultP);
	INFO("N+1 sector: ||krylovExpmvCSR(e_0)|| = " << normP);
	CHECK(std::abs(normP - 1.0) < 1e-10);
}

// ── TC-7: Gα/Gβ spin-seed averaging (GBEK Eq. 70) ──────────────────────────────
//
// Regression test for a real bug (found and fixed 2026-07-08): a
// spin-imbalanced base filling (nup != ndown, e.g. the true atomic-limit
// single-atom seed nup=1, ndown=0) produces a spin-polarized extended-Fock
// ground state on its own. Without averaging Galpha (impurity's extra
// electron up) and Gbeta (down), the up-channel occupation is frozen at
// n=1 for all time instead of oscillating around the physically-required
// n=1/2 -- exactly the failure this test guards against. The expected
// value 0.5 is not a tolerance-tuned regression anchor: it follows directly
// from GBEK Sec. VI-A (paramagnetic, particle-hole symmetric state at half
// filling, <n_sigma>=1/2 for all t) and from the exact self-consistency
// relation Delta(t,t)=v(t)^2 n(t) matching the paper's own Fig. 3 colorbar
// maximum of 0.5.
//
// nBath=0 (true atomic limit, Delta^-=0 exactly) isolates this from the
// separate second-bath Cholesky machinery: with an empty bathParams vector,
// solveLplus builds nsites_ext=1+2L with no first-bath sites at all, so any
// deviation from n=1/2 can only come from the (nup_ext,ndown_ext) seeding,
// not from Cholesky-approximation error.
TEST_CASE("Gimp(t,t) plateaus at 1/2 for a spin-imbalanced atomic-limit seed",
          "[GBEK][spin-averaging]")
{
	static const std::string kAtomicConfig
	    = "##Ainur1.0\n\n"
	      "FicticiousBeta=10;\n"
	      "ChemicalPotential=0.;\n"
	      "Matsubaras=20;\n"
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
	      "TargetElectronsDown=0;\n"
	      "int ImpuritySite=0;\n"
	      "real HubbardU=2.;\n"
	      "HubbardUFinal=2.;\n"
	      "RootOutputname=\"testGBEKAtomic\";\n"
	      "InfiniteLoopKeptStates=20;\n"
	      "matrix FiniteLoopsGs=[[@auto, 20, 0],[@auto, 20, 0]];\n"
	      "real OmegaBegin=-4.;\n"
	      "integer OmegaTotal=40;\n"
	      "real OmegaStep=0.2;\n"
	      "real OmegaDelta=0.2;\n"
	      "integer TridiagSteps=100;\n"
	      "real TridiagEps=1e-6;\n"
	      "TruncationTolerance=\"1e-6,20\";\n"
	      "CorrectionVectorEta=0.;\n"
	      "GsWeight=0.1;\n"
	      "matrix FiniteLoopsOmega=[[@auto, 20, 2],[@auto, 20, 2]];\n"
	      "TmaxNeq=0.2;\n"
	      "NtNeq=4;\n"
	      "NeqDmftIter=1;\n"
	      "NeqDmftTolerance=1e-4;\n"
	      "NeqSolver=\"gbek\";\n"
	      "NeqBathRank=2;\n"
	      "BandwidthFinal=4.;\n";

	InputNgType::Writeable ioW(Dmft::CincuentaInputCheck {}, kAtomicConfig);
	InputNgType::Readable  io(ioW);
	ParamsType             params(io);
	SolverType             solver(params, io);

	solver.solve({}); // empty bathParams: nBath=0, true atomic limit

	using KBType = SolverType::KBType;
	KBType slice(params.nT,
	             params.eqParams.nMatsubaras,
	             params.dt,
	             params.eqParams.ficticiousBeta
	                 / static_cast<RealType>(params.eqParams.nMatsubaras));

	for (int n = 0; n <= static_cast<int>(params.nT); ++n) {
		solver.computeGimp(slice, n);
		const ComplexType g
		    = slice.lesser(static_cast<SizeType>(n), static_cast<SizeType>(n));
		INFO("n=" << n << " G_imp(t,t)=" << g);
		CHECK(g.real() == Catch::Approx(0.0).margin(1e-10));
		CHECK(g.imag() == Catch::Approx(0.5).epsilon(1e-8));
	}
}

// ── Hamiltonian-construction regression check ────────────────────────────────
// Added 2026-07-13 while investigating a "wrong phase in the reconstructed
// Weiss field" bug (see cincuenta/TestSuite/gbek_reference/README.md). Row-
// by-row comparisons of the Cholesky factor V, or end-to-end comparisons of
// the converged Weiss field, cannot cleanly localize a Hamiltonian bug: V
// and the converged Weiss field both depend on the whole self-consistency
// loop's history and are not gauge/basis invariant, so a mismatch there
// could come from many places. Eigenvalues of a single, FIXED (externally
// supplied, not self-consistently derived) V are different: they are
// basis-ordering-independent, so comparing them against an INDEPENDENT
// from-scratch reconstruction (cross_check_gbek_hamiltonian.py in
// gbek_reference, which never calls into cincuenta or LanczosPlusPlus) can
// only disagree because of the Hamiltonian's actual physics -- wrong
// Jordan-Wigner sign, wrong site layout, wrong occ/empty coupling polarity,
// or a wrong diagonal/potential term. This test hardcodes that Python
// script's output; rerun the script (see its own docstring) if the
// Hamiltonian construction is ever intentionally changed, and update the
// expected values here to match.
TEST_CASE("Hamiltonian eigenvalue spectrum matches an independent Python "
          "reconstruction (nBath=0, L=1, complex vMid)",
          "[GBEK][Hamiltonian]")
{
	// Same small system as the "Gimp(t,t) plateaus..." test just above:
	// nBath=0 (true atomic limit, passed as an empty bathParams vector),
	// L=1 second-bath pair, spin-imbalanced nup=1/ndown=0 seed.
	static const std::string kConfig
	    = "##Ainur1.0\n\n"
	      "FicticiousBeta=10;\n"
	      "ChemicalPotential=0.;\n"
	      "Matsubaras=20;\n"
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
	      "TargetElectronsDown=0;\n"
	      "int ImpuritySite=0;\n"
	      "real HubbardU=2.;\n"
	      "HubbardUFinal=2.;\n"
	      "RootOutputname=\"testGBEKHamiltonian\";\n"
	      "InfiniteLoopKeptStates=20;\n"
	      "matrix FiniteLoopsGs=[[@auto, 20, 0],[@auto, 20, 0]];\n"
	      "real OmegaBegin=-4.;\n"
	      "integer OmegaTotal=40;\n"
	      "real OmegaStep=0.2;\n"
	      "real OmegaDelta=0.2;\n"
	      "integer TridiagSteps=100;\n"
	      "real TridiagEps=1e-6;\n"
	      "TruncationTolerance=\"1e-6,20\";\n"
	      "CorrectionVectorEta=0.;\n"
	      "GsWeight=0.1;\n"
	      "matrix FiniteLoopsOmega=[[@auto, 20, 2],[@auto, 20, 2]];\n"
	      "TmaxNeq=0.2;\n"
	      "NtNeq=2;\n"
	      "NeqDmftIter=1;\n"
	      "NeqDmftTolerance=1e-4;\n"
	      "NeqSolver=\"gbek\";\n"
	      "NeqBathRank=1;\n"
	      "BandwidthFinal=4.;\n";

	InputNgType::Writeable ioW(Dmft::CincuentaInputCheck {}, kConfig);
	InputNgType::Readable  io(ioW);
	ParamsType             params(io);
	SolverType             solver(params, io);
	solver.solve({}); // empty bathParams: nBath=0

	// Matches cross_check_gbek_hamiltonian.py's vMid exactly.
	const std::vector<ComplexType> vMid = { ComplexType(0.3, 0.2) };

	auto Hnm1 = buildDenseFromApplyHext(
	    solver, vMid, Acc::upWordsNm1(solver), Acc::dnWordsNm1(solver), Acc::dim1Nm1(solver));
	REQUIRE(Hnm1.n_row() == 9);
	assertHermitian(Hnm1);

	PsimagLite::Vector<RealType>::Type eigs;
	PsimagLite::diag(Hnm1, eigs, 'N');
	std::sort(eigs.begin(), eigs.end());

	// From cross_check_gbek_hamiltonian.py (independent, from-scratch
	// reconstruction -- see that script's docstring).
	const std::vector<RealType> expected
	    = { -1.63578166916006, -1.21414284285429, -1.21414284285429, -1.0, 0.0, 0.0,
		0.21414284285429,  0.21414284285429,  0.63578166916006 };
	REQUIRE(eigs.size() == expected.size());
	for (SizeType i = 0; i < expected.size(); ++i) {
		INFO("eigenvalue index " << i << ": got " << eigs[i] << ", expected "
		                         << expected[i]);
		CHECK(eigs[i] == Catch::Approx(expected[i]).margin(1e-10));
	}
}

// ── Propagation regression check ─────────────────────────────────────────────
// Added 2026-07-13, alongside the Hamiltonian eigenvalue check above. That
// test only exercises a single, static Hamiltonian for one fixed V; it says
// nothing about the seeded initial state (seedState) or the Krylov time
// propagation itself. This test covers both: it propagates the SAME small
// system two steps with two DIFFERENT, externally-fixed vMid values
// (bypassing NeqBathDecomposition and the self-consistency loop entirely,
// exactly like the Hamiltonian check), and compares the resulting G^<(n,j)
// against an independent from-scratch reconstruction
// (cross_check_gbek_propagation.py in gbek_reference) that starts from the
// SAME seeded initial state -- confirmed by direct inspection to be the
// pure basis state with both second-bath sites in their t=0 configuration
// and the impurity/empty site empty (see that script's docstring for the
// derivation). Rerun the Python script and update the expected values here
// if the propagation or seeding is ever intentionally changed.
TEST_CASE("Time propagation matches an independent Python reconstruction "
          "(nBath=0, L=1, two-step complex vMid)",
          "[GBEK][Krylov][propagation]")
{
	static const std::string kConfig
	    = "##Ainur1.0\n\n"
	      "FicticiousBeta=10;\nChemicalPotential=0.;\nMatsubaras=20;\n"
	      "LatticeGf=\"energy,semicircular,4\";\nNumberOfBathPoints=1;\n"
	      "DmftNumberOfIterations=1;\nDmftTolerance=1e-3;\n"
	      "ImpuritySolver=\"exactdiag\";\nFitOptions=particleholesymmetric;\n"
	      "MinParamsDelta=0.01;\nMinParamsDelta2=0.01;\nMinParamsTolerance=1e-4;\n"
	      "MinParamsMaxIter=100;\nMinParamsVerbose=0;\n"
	      "TargetElectronsUp=1;\nTargetElectronsDown=0;\n"
	      "int ImpuritySite=0;\nreal HubbardU=2.;\nHubbardUFinal=2.;\n"
	      "RootOutputname=\"testGBEKPropagation\";\nInfiniteLoopKeptStates=20;\n"
	      "matrix FiniteLoopsGs=[[@auto, 20, 0],[@auto, 20, 0]];\n"
	      "real OmegaBegin=-4.;\ninteger OmegaTotal=40;\nreal OmegaStep=0.2;\n"
	      "real OmegaDelta=0.2;\ninteger TridiagSteps=100;\nreal TridiagEps=1e-6;\n"
	      "TruncationTolerance=\"1e-6,20\";\nCorrectionVectorEta=0.;\nGsWeight=0.1;\n"
	      "matrix FiniteLoopsOmega=[[@auto, 20, 2],[@auto, 20, 2]];\n"
	      "TmaxNeq=0.2;\nNtNeq=2;\nNeqDmftIter=1;\nNeqDmftTolerance=1e-4;\n"
	      "NeqSolver=\"gbek\";\nNeqBathRank=1;\nBandwidthFinal=4.;\n";

	InputNgType::Writeable ioW(Dmft::CincuentaInputCheck {}, kConfig);
	InputNgType::Readable  io(ioW);
	ParamsType             params(io);
	SolverType             solver(params, io);
	solver.solve({}); // empty bathParams: nBath=0

	// Matches cross_check_gbek_propagation.py's dt and two vMid values exactly.
	const RealType    dt = 0.05;
	const ComplexType vMidStep0(0.3, 0.2);
	const ComplexType vMidStep1(0.1, -0.15);

	CrsMatrixType& csr = Acc::csrNm1Mut(solver);

	Acc::updateCSR(solver, csr, Acc::varNm1(solver), { vMidStep0 });
	VectorComplex psi1 = Acc::krylovExpmvCSR(solver, Acc::bStates(solver)[0], csr, dt);

	Acc::updateCSR(solver, csr, Acc::varNm1(solver), { vMidStep1 });
	VectorComplex psi2 = Acc::krylovExpmvCSR(solver, psi1, csr, dt);

	const std::vector<VectorComplex> psis = { Acc::bStates(solver)[0], psi1, psi2 };

	auto innerProd = [](const VectorComplex& a, const VectorComplex& b)
	{
		ComplexType s(0);
		for (SizeType i = 0; i < a.size(); ++i)
			s += std::conj(a[i]) * b[i];
		return s;
	};

	for (SizeType k = 0; k < psis.size(); ++k) {
		INFO("norm at step " << k);
		CHECK(std::real(innerProd(psis[k], psis[k])) == Catch::Approx(1.0).margin(1e-10));
	}

	// From cross_check_gbek_propagation.py: G_lesser(n,j) = i * conj(psi_j).psi_n,
	// for (n,j) with j <= n.
	struct Expected {
		int         n, j;
		ComplexType value;
	};
	const std::vector<Expected> expected = {
		{ 0, 0, ComplexType(0.00000000000000, 1.00000000000000) },
		{ 1, 0, ComplexType(0.00000541510959, 0.99967511169887) },
		{ 1, 1, ComplexType(0.00000000000000, 1.00000000000000) },
		{ 2, 0, ComplexType(0.00033118585217, 0.99961010253565) },
		{ 2, 1, ComplexType(0.00061824086200, 0.99993403666773) },
		{ 2, 2, ComplexType(0.00000000000000, 1.00000000000000) },
	};
	for (const auto& e : expected) {
		const ComplexType g = ComplexType(0, 1) * innerProd(psis[e.j], psis[e.n]);
		INFO("G_lesser(n=" << e.n << ",j=" << e.j << "): got " << g << ", expected "
		                   << e.value);
		CHECK(g.real() == Catch::Approx(e.value.real()).margin(1e-8));
		CHECK(g.imag() == Catch::Approx(e.value.imag()).margin(1e-8));
	}
}

// ── Seed-scheme regression check ─────────────────────────────────────────────
// Added 2026-07-13, guarding against a real bug found and fixed the same day:
// ImpuritySolverNeqGBEK.h used to seed c_{imp,up} ONCE (at t=0, from the
// pre-quench GS) and propagate the resulting (N-1)-sector state forward using
// only the (N-1)-sector Hamiltonian. That is mathematically wrong whenever
// c does not commute with H (true here: [U n_up n_dn, c_up] != 0) -- the
// correct Heisenberg-picture construction requires re-seeding
// c_{imp,up}/c†_{imp,up} against the propagated N-sector reference state
// PhiNHist AT EVERY time step, not just at t=0 (see propagateOneStep's and
// gLesserRowGBEKSector's doc comments in ImpuritySolverNeqGBEK.h). The old,
// wrong scheme produced a systematically suppressed off-diagonal (two-time)
// imaginary part, growing with |n-j|, that this test would have caught.
//
// This test drives the production N-sector CSR (csrN) and the c_{imp,up}
// operator (cUpNm1) directly -- bypassing NeqBathDecomposition/the
// self-consistency loop entirely, exactly like the propagation test above --
// with a real, time-dependent cosine-ramp hopping (not an arbitrary fixed
// complex vMid), so this exercises the ACTUAL reseed-at-every-step code path
// used in production, not a hand reimplementation of it.
//
// Expected values are G_bruteforce from
// cincuenta/TestSuite/gbek_reference/cross_check_seed_scheme.py, an
// independent, from-scratch Python reconstruction using dense
// scipy.linalg.expm (no shared code with gbek_dynamics.py's propagate() or
// compute_g_lesser(), so it cannot inherit a shared bug). That script's
// docstring documents the full derivation and how its psi0 was verified to
// match cincuenta's actual pre-quench seed; rerun it and update the expected
// values here if the seeding or propagation is ever intentionally changed.
TEST_CASE("G^< two-time construction matches an independent Python "
          "brute-force ground truth (nBath=0, L=1, reseed-at-every-step)",
          "[GBEK][seed-scheme]")
{
	static const std::string kConfig
	    = "##Ainur1.0\n\n"
	      "FicticiousBeta=10;\nChemicalPotential=0.;\nMatsubaras=20;\n"
	      "LatticeGf=\"energy,semicircular,4\";\nNumberOfBathPoints=1;\n"
	      "DmftNumberOfIterations=1;\nDmftTolerance=1e-3;\n"
	      "ImpuritySolver=\"exactdiag\";\nFitOptions=particleholesymmetric;\n"
	      "MinParamsDelta=0.01;\nMinParamsDelta2=0.01;\nMinParamsTolerance=1e-4;\n"
	      "MinParamsMaxIter=100;\nMinParamsVerbose=0;\n"
	      "TargetElectronsUp=1;\nTargetElectronsDown=0;\n"
	      "int ImpuritySite=0;\nreal HubbardU=2.;\nHubbardUFinal=2.;\n"
	      "RootOutputname=\"testGBEKSeedScheme\";\nInfiniteLoopKeptStates=20;\n"
	      "matrix FiniteLoopsGs=[[@auto, 20, 0],[@auto, 20, 0]];\n"
	      "real OmegaBegin=-4.;\ninteger OmegaTotal=40;\nreal OmegaStep=0.2;\n"
	      "real OmegaDelta=0.2;\ninteger TridiagSteps=100;\nreal TridiagEps=1e-6;\n"
	      "TruncationTolerance=\"1e-6,20\";\nCorrectionVectorEta=0.;\nGsWeight=0.1;\n"
	      "matrix FiniteLoopsOmega=[[@auto, 20, 2],[@auto, 20, 2]];\n"
	      "TmaxNeq=0.3;\nNtNeq=6;\nNeqDmftIter=1;\nNeqDmftTolerance=1e-4;\n"
	      "NeqSolver=\"gbek\";\nNeqBathRank=1;\nBandwidthFinal=4.;\n";

	InputNgType::Writeable ioW(Dmft::CincuentaInputCheck {}, kConfig);
	InputNgType::Readable  io(ioW);
	ParamsType             params(io);
	SolverType             solver(params, io);
	solver.solve({}); // empty bathParams: nBath=0

	// Matches cross_check_seed_scheme.py's dt, tq, tstar_f, nsteps exactly.
	const RealType dt      = 0.05;
	const RealType tq      = 0.25;
	const RealType tstar_f = 1.0;
	const int      nsteps  = 6;
	const RealType pi      = RealType(3.14159265358979323846);

	auto vRamp = [&](RealType t) -> RealType
	{
		if (t <= 0)
			return 0.0;
		if (t >= tq)
			return 1.0;
		return 0.5 * (1.0 - std::cos(pi * t / tq));
	};
	auto hop = [&](RealType t) { return tstar_f * vRamp(t); };

	// Forward-propagate the N-sector reference trajectory PhiNHist, then
	// reseed c_{imp,up} at every step to get bStates[n] -- the ACTUAL
	// production construction (ImpuritySolverNeqGBEK.h::propagateOneStep),
	// driven here with an externally-fixed hop(t) schedule instead of
	// NeqBathDecomposition's self-consistent V.
	std::vector<VectorComplex> phiN(static_cast<SizeType>(nsteps + 1));
	phiN[0] = Acc::phiNHist(solver)[0];

	CrsMatrixType& csrN = Acc::csrNMut(solver);
	for (int n = 1; n <= nsteps; ++n) {
		const RealType tMid = (n - 0.5) * dt;
		Acc::updateCSR(solver, csrN, Acc::varN(solver), { ComplexType(hop(tMid)) });
		phiN[static_cast<SizeType>(n)]
		    = Acc::krylovExpmvCSR(solver, phiN[static_cast<SizeType>(n - 1)], csrN, dt);
	}

	std::vector<VectorComplex> bStatesLocal(static_cast<SizeType>(nsteps + 1));
	for (int n = 0; n <= nsteps; ++n)
		bStatesLocal[static_cast<SizeType>(n)] = Acc::sparseMatVec(
		    solver, Acc::cUpNm1(solver), phiN[static_cast<SizeType>(n)]);

	// Sanity check: bStatesLocal[0] must match the production seed exactly
	// (both are c_{imp,up}|pre-quench GS>, just reached via this test's
	// explicit sparseMatVec vs. the solver's own internal seeding).
	{
		const auto& prodSeed = Acc::bStates(solver)[0];
		REQUIRE(prodSeed.size() == bStatesLocal[0].size());
		RealType maxDiff = 0;
		for (SizeType i = 0; i < prodSeed.size(); ++i)
			maxDiff = std::max(maxDiff, std::abs(prodSeed[i] - bStatesLocal[0][i]));
		INFO("max|bStates[0] - locally-derived seed| = " << maxDiff);
		CHECK(maxDiff < 1e-12);
	}

	// G^<(n,j) row via backward Krylov sweep from bStatesLocal[n] -- mirrors
	// gLesserRowGBEKSector exactly, driven with the same externally-fixed
	// hop(t) schedule.
	auto innerProd = [](const VectorComplex& a, const VectorComplex& b)
	{
		ComplexType s(0);
		for (SizeType i = 0; i < a.size(); ++i)
			s += std::conj(a[i]) * b[i];
		return s;
	};

	CrsMatrixType& csrNm1     = Acc::csrNm1Mut(solver);
	auto           gLesserRow = [&](int n)
	{
		std::vector<ComplexType> row(static_cast<SizeType>(n + 1));
		VectorComplex            psi  = bStatesLocal[static_cast<SizeType>(n)];
		row[static_cast<SizeType>(n)] = ComplexType(0, 1) * innerProd(psi, psi);
		for (int k = n - 1; k >= 0; --k) {
			const RealType tMid = (k + 0.5) * dt;
			Acc::updateCSR(
			    solver, csrNm1, Acc::varNm1(solver), { ComplexType(hop(tMid)) });
			psi = Acc::krylovExpmvCSR(solver, psi, csrNm1, -dt);
			row[static_cast<SizeType>(k)] = ComplexType(0, 1)
			    * innerProd(bStatesLocal[static_cast<SizeType>(k)], psi);
		}
		return row;
	};

	// From cross_check_seed_scheme.py's G_bruteforce array (see that
	// script's docstring for the full derivation).
	struct Expected {
		int         n, j;
		ComplexType value;
	};
	const std::vector<Expected> expected = {
		{ 1, 0, ComplexType(-0.049979, 0.998749) },
		{ 1, 1, ComplexType(0.000000, 0.999999) },
		{ 2, 0, ComplexType(-0.099818, 0.994872) },
		{ 2, 1, ComplexType(-0.049971, 0.998630) },
		{ 2, 2, ComplexType(0.000000, 0.999867) },
		{ 3, 0, ComplexType(-0.149216, 0.987452) },
		{ 3, 1, ComplexType(-0.099678, 0.993719) },
		{ 3, 2, ComplexType(-0.049890, 0.997705) },
		{ 3, 3, ComplexType(0.000000, 0.998668) },
		{ 4, 0, ComplexType(-0.197386, 0.974376) },
		{ 4, 1, ComplexType(-0.148443, 0.983115) },
		{ 4, 2, ComplexType(-0.099136, 0.989964) },
		{ 4, 3, ComplexType(-0.049557, 0.994403) },
		{ 4, 4, ComplexType(0.000000, 0.994222) },
		{ 5, 0, ComplexType(-0.243052, 0.953805) },
		{ 5, 1, ComplexType(-0.195085, 0.964912) },
		{ 5, 2, ComplexType(-0.146661, 0.974618) },
		{ 5, 3, ComplexType(-0.097823, 0.982690) },
		{ 5, 4, ComplexType(-0.048718, 0.986937) },
		{ 5, 5, ComplexType(0.000000, 0.984579) },
		{ 6, 0, ComplexType(-0.285118, 0.926218) },
		{ 6, 1, ComplexType(-0.238484, 0.939522) },
		{ 6, 2, ComplexType(-0.191328, 0.951939) },
		{ 6, 3, ComplexType(-0.143657, 0.963540) },
		{ 6, 4, ComplexType(-0.095472, 0.972174) },
		{ 6, 5, ComplexType(-0.047250, 0.974792) },
		{ 6, 6, ComplexType(0.000000, 0.970117) },
	};

	std::vector<std::vector<ComplexType>> rows(static_cast<SizeType>(nsteps + 1));
	for (int n = 1; n <= nsteps; ++n)
		rows[static_cast<SizeType>(n)] = gLesserRow(n);

	for (const auto& e : expected) {
		const ComplexType g = rows[static_cast<SizeType>(e.n)][static_cast<SizeType>(e.j)];
		INFO("G_lesser(n=" << e.n << ",j=" << e.j << "): got " << g << ", expected "
		                   << e.value);
		CHECK(g.real() == Catch::Approx(e.value.real()).margin(1e-5));
		CHECK(g.imag() == Catch::Approx(e.value.imag()).margin(1e-5));
	}
}

// ── Double occupation / energy observable tests (paper Figs. 9-10) ──────────
// Added 2026-07-15, alongside the C++ port of the Python reference's Fig.
// 9/10 observables (see project memory project_gbek_ekin_factor2_bugfix and
// project_gbek_double_occupation_next). Both tests below hand-drive PhiNHist
// directly (same style as the seed-scheme test above), bypassing
// NeqBathDecomposition/self-consistency entirely, so they are fast and
// deterministic.

// ── TC: trivial zero-coupling anchor ─────────────────────────────────────────
// V(t)==0 for all t: the impurity/bath coupling is never switched on, so
// PhiNHist never leaves its initial (diagonal-Hamiltonian) eigenstate up to
// an overall phase. docc(t) and Ekin(t) must be EXACTLY 0 for all t, and
// Etot(t) exactly -U/4 -- this is a basic plumbing/indexing sanity check
// (diagFixed captured correctly, docc's basis-index decomposition correct),
// NOT a regression guard for the factor-of-2 bug: with V==0, <H_hop>==0
// identically regardless of any overall scale factor applied to it, so a
// 2x (or any other multiplicative) bug would be invisible here. See the
// conservation-law test below for the check that actually catches that bug.
TEST_CASE("docc/Ekin are exactly zero when the second-bath coupling is never "
          "switched on",
          "[GBEK][energy][docc]")
{
	static const std::string kConfig
	    = "##Ainur1.0\n\n"
	      "FicticiousBeta=10;\nChemicalPotential=0.;\nMatsubaras=20;\n"
	      "LatticeGf=\"energy,semicircular,4\";\nNumberOfBathPoints=1;\n"
	      "DmftNumberOfIterations=1;\nDmftTolerance=1e-3;\n"
	      "ImpuritySolver=\"exactdiag\";\nFitOptions=particleholesymmetric;\n"
	      "MinParamsDelta=0.01;\nMinParamsDelta2=0.01;\nMinParamsTolerance=1e-4;\n"
	      "MinParamsMaxIter=100;\nMinParamsVerbose=0;\n"
	      "TargetElectronsUp=1;\nTargetElectronsDown=0;\n"
	      "int ImpuritySite=0;\nreal HubbardU=2.;\nHubbardUFinal=2.;\n"
	      "RootOutputname=\"testGBEKZeroCoupling\";\nInfiniteLoopKeptStates=20;\n"
	      "matrix FiniteLoopsGs=[[@auto, 20, 0],[@auto, 20, 0]];\n"
	      "real OmegaBegin=-4.;\ninteger OmegaTotal=40;\nreal OmegaStep=0.2;\n"
	      "real OmegaDelta=0.2;\ninteger TridiagSteps=100;\nreal TridiagEps=1e-6;\n"
	      "TruncationTolerance=\"1e-6,20\";\nCorrectionVectorEta=0.;\nGsWeight=0.1;\n"
	      "matrix FiniteLoopsOmega=[[@auto, 20, 2],[@auto, 20, 2]];\n"
	      "TmaxNeq=0.5;\nNtNeq=10;\nNeqDmftIter=1;\nNeqDmftTolerance=1e-4;\n"
	      "NeqSolver=\"gbek\";\nNeqBathRank=1;\nBandwidthFinal=4.;\n";

	InputNgType::Writeable ioW(Dmft::CincuentaInputCheck {}, kConfig);
	InputNgType::Readable  io(ioW);
	ParamsType             params(io);
	SolverType             solver(params, io);
	solver.solve({}); // empty bathParams: nBath=0

	const RealType dt     = 0.05;
	const int      nsteps = 10;
	const RealType uFinal = 2.0;

	std::vector<VectorComplex>& phiN = Acc::phiNHistMut(solver);
	CrsMatrixType&              csrN = Acc::csrNMut(solver);
	for (int n = 1; n <= nsteps; ++n) {
		Acc::updateCSR(solver, csrN, Acc::varN(solver), { ComplexType(0, 0) });
		phiN[static_cast<SizeType>(n)]
		    = Acc::krylovExpmvCSR(solver, phiN[static_cast<SizeType>(n - 1)], csrN, dt);
	}

	for (int n = 0; n <= nsteps; ++n) {
		const RealType docc = Acc::computeDoubleOccupationGBEKSector(solver, n);
		const RealType ekin = Acc::computeKineticEnergyGBEKSector(solver, n);
		const RealType eint = uFinal * (docc - 0.25);
		const RealType etot = ekin + eint;
		INFO("n=" << n << " docc=" << docc << " ekin=" << ekin << " etot=" << etot);
		CHECK(docc == Catch::Approx(0.0).margin(1e-12));
		CHECK(ekin == Catch::Approx(0.0).margin(1e-10));
		CHECK(etot == Catch::Approx(-uFinal / 4).margin(1e-10));
	}
}

// ── TC: energy-transfer plausibility check (NOT a factor-of-2 regression
// guard -- see the dense cross-check test below for that) ───────────────────
// Drives PhiNHist with the SAME physical cosine ramp (tq=0.25, tstar_f=1.0)
// used by the actual GBEK Fig. 9 protocol (identical to the seed-scheme
// test's hop(t) above), at nBath=0/L=1, and checks Etot(t) stays roughly
// close to its analytic t=0 value -U/4 shortly after the ramp ("tiny energy
// transfer", paper Sec. VI-C).
//
// IMPORTANT LIMITATION, found empirically while writing this test: once
// hop(t) becomes constant (t > t_q), <psi(t)|H(t_q)|psi(t)> is EXACTLY
// conserved by plain unitarity, for ANY multiplicative scale applied to how
// that number gets reported as "Ekin" -- a 2x-too-large Ekin is just as
// "conserved" as the correct one, only at a different (still plausible)
// plateau value. Verified directly: temporarily reintroducing the
// factor-of-2 bug in ImpuritySolverNeqGBEK.h and rerunning this exact test
// produced Etot values of -0.506 (bug) vs. -0.474..-0.393 (fixed) at
// n=6..10 -- BOTH within the tolerance below, i.e. this test would NOT have
// caught the bug. It is kept as a basic plausibility check (values are
// finite and roughly the right order of magnitude), not a regression guard.
TEST_CASE("Etot(t) stays roughly close to its analytic -U/4 value shortly "
          "after the ramp completes (plausibility check only)",
          "[GBEK][energy][conservation]")
{
	static const std::string kConfig
	    = "##Ainur1.0\n\n"
	      "FicticiousBeta=10;\nChemicalPotential=0.;\nMatsubaras=20;\n"
	      "LatticeGf=\"energy,semicircular,4\";\nNumberOfBathPoints=1;\n"
	      "DmftNumberOfIterations=1;\nDmftTolerance=1e-3;\n"
	      "ImpuritySolver=\"exactdiag\";\nFitOptions=particleholesymmetric;\n"
	      "MinParamsDelta=0.01;\nMinParamsDelta2=0.01;\nMinParamsTolerance=1e-4;\n"
	      "MinParamsMaxIter=100;\nMinParamsVerbose=0;\n"
	      "TargetElectronsUp=1;\nTargetElectronsDown=0;\n"
	      "int ImpuritySite=0;\nreal HubbardU=2.;\nHubbardUFinal=2.;\n"
	      "RootOutputname=\"testGBEKConservation\";\nInfiniteLoopKeptStates=20;\n"
	      "matrix FiniteLoopsGs=[[@auto, 20, 0],[@auto, 20, 0]];\n"
	      "real OmegaBegin=-4.;\ninteger OmegaTotal=40;\nreal OmegaStep=0.2;\n"
	      "real OmegaDelta=0.2;\ninteger TridiagSteps=100;\nreal TridiagEps=1e-6;\n"
	      "TruncationTolerance=\"1e-6,20\";\nCorrectionVectorEta=0.;\nGsWeight=0.1;\n"
	      "matrix FiniteLoopsOmega=[[@auto, 20, 2],[@auto, 20, 2]];\n"
	      "TmaxNeq=1.5;\nNtNeq=30;\nNeqDmftIter=1;\nNeqDmftTolerance=1e-4;\n"
	      "NeqSolver=\"gbek\";\nNeqBathRank=1;\nBandwidthFinal=4.;\n";

	InputNgType::Writeable ioW(Dmft::CincuentaInputCheck {}, kConfig);
	InputNgType::Readable  io(ioW);
	ParamsType             params(io);
	SolverType             solver(params, io);
	solver.solve({}); // empty bathParams: nBath=0

	const RealType dt      = 0.05;
	const RealType tq      = 0.25;
	const RealType tstar_f = 1.0;
	const int      nsteps  = 30;
	const RealType uFinal  = 2.0;
	const RealType pi      = RealType(3.14159265358979323846);

	auto vRamp = [&](RealType t) -> RealType
	{
		if (t <= 0)
			return 0.0;
		if (t >= tq)
			return 1.0;
		return 0.5 * (1.0 - std::cos(pi * t / tq));
	};
	auto hop = [&](RealType t) { return tstar_f * vRamp(t); };

	std::vector<VectorComplex>& phiN = Acc::phiNHistMut(solver);
	CrsMatrixType&              csrN = Acc::csrNMut(solver);
	for (int n = 1; n <= nsteps; ++n) {
		const RealType tMid = (n - 0.5) * dt;
		Acc::updateCSR(solver, csrN, Acc::varN(solver), { ComplexType(hop(tMid)) });
		phiN[static_cast<SizeType>(n)]
		    = Acc::krylovExpmvCSR(solver, phiN[static_cast<SizeType>(n - 1)], csrN, dt);
	}

	std::vector<RealType> etot(static_cast<SizeType>(nsteps + 1));
	for (int n = 0; n <= nsteps; ++n) {
		const RealType docc            = Acc::computeDoubleOccupationGBEKSector(solver, n);
		const RealType ekin            = Acc::computeKineticEnergyGBEKSector(solver, n);
		const RealType eint            = uFinal * (docc - 0.25);
		etot[static_cast<SizeType>(n)] = ekin + eint;
	}

	// n=0: exact analytic value, no approximation involved.
	CHECK(etot[0] == Catch::Approx(-uFinal / 4).margin(1e-10));

	// Shortly after the ramp completes (tq=0.25 -> n=5), Etot should still
	// be close to -U/4 ("tiny energy transfer"). Loose margin: this is an
	// L=1 toy instance of the real physics, not the converged Lbath=8 case.
	const RealType tol = 0.3;
	for (int n = 6; n <= 10; ++n) {
		INFO("n=" << n << " t=" << n * dt << " Etot=" << etot[static_cast<SizeType>(n)]);
		CHECK(etot[static_cast<SizeType>(n)] == Catch::Approx(-uFinal / 4).margin(tol));
	}
}

// ── TC: dense cross-check (the ACTUAL regression guard for the factor-of-2
// bug) ────────────────────────────────────────────────────────────────────
// Independently recomputes <Ekin> for a single hand-picked step via the
// DENSE applyHext path (buildDenseFromApplyHext, already cross-validated
// against the production CSR Hamiltonian by "buildHextCSR produces
// Hermitian matrix..." and "CSR SpMV agrees with applyHext..." above), with
// the 0.5 factor and onsite-diagonal subtraction written out explicitly
// HERE rather than delegated to computeKineticEnergyGBEKSector -- so a
// future change to that function's 0.5 factor (or its diagFixed usage) that
// isn't ALSO made here will show up as a mismatch, not a silent pass.
TEST_CASE("computeKineticEnergyGBEKSector matches an independent dense-matrix "
          "computation of 0.5*(<H_hop>-onsite)",
          "[GBEK][energy][dense-cross-check]")
{
	static const std::string kConfig
	    = "##Ainur1.0\n\n"
	      "FicticiousBeta=10;\nChemicalPotential=0.;\nMatsubaras=20;\n"
	      "LatticeGf=\"energy,semicircular,4\";\nNumberOfBathPoints=1;\n"
	      "DmftNumberOfIterations=1;\nDmftTolerance=1e-3;\n"
	      "ImpuritySolver=\"exactdiag\";\nFitOptions=particleholesymmetric;\n"
	      "MinParamsDelta=0.01;\nMinParamsDelta2=0.01;\nMinParamsTolerance=1e-4;\n"
	      "MinParamsMaxIter=100;\nMinParamsVerbose=0;\n"
	      "TargetElectronsUp=1;\nTargetElectronsDown=0;\n"
	      "int ImpuritySite=0;\nreal HubbardU=2.;\nHubbardUFinal=2.;\n"
	      "RootOutputname=\"testGBEKDenseCrossCheck\";\nInfiniteLoopKeptStates=20;\n"
	      "matrix FiniteLoopsGs=[[@auto, 20, 0],[@auto, 20, 0]];\n"
	      "real OmegaBegin=-4.;\ninteger OmegaTotal=40;\nreal OmegaStep=0.2;\n"
	      "real OmegaDelta=0.2;\ninteger TridiagSteps=100;\nreal TridiagEps=1e-6;\n"
	      "TruncationTolerance=\"1e-6,20\";\nCorrectionVectorEta=0.;\nGsWeight=0.1;\n"
	      "matrix FiniteLoopsOmega=[[@auto, 20, 2],[@auto, 20, 2]];\n"
	      "TmaxNeq=0.2;\nNtNeq=2;\nNeqDmftIter=1;\nNeqDmftTolerance=1e-4;\n"
	      "NeqSolver=\"gbek\";\nNeqBathRank=1;\nBandwidthFinal=4.;\n";

	InputNgType::Writeable ioW(Dmft::CincuentaInputCheck {}, kConfig);
	InputNgType::Readable  io(ioW);
	ParamsType             params(io);
	SolverType             solver(params, io);
	solver.solve({}); // empty bathParams: nBath=0

	const RealType    dt = 0.05;
	const ComplexType vMid(0.3, 0.2); // arbitrary, matches other tests' convention

	std::vector<VectorComplex>& phiN = Acc::phiNHistMut(solver);
	CrsMatrixType&              csrN = Acc::csrNMut(solver);
	Acc::updateCSR(solver, csrN, Acc::varN(solver), { vMid });
	phiN[1] = Acc::krylovExpmvCSR(solver, phiN[0], csrN, dt);

	// Independent dense H(vMid) for the N sector -- a different code path
	// (applyHext, not the CSR csrN this test just drove) from the same
	// underlying physics, already cross-validated elsewhere in this file.
	auto denseH = buildDenseFromApplyHext(
	    solver, { vMid }, Acc::upWordsN(solver), Acc::dnWordsN(solver), Acc::dim1N(solver));

	const VectorComplex& psi1 = phiN[1];
	const SizeType       dim  = psi1.size();
	ComplexType          raw(0);
	for (SizeType i = 0; i < dim; ++i) {
		ComplexType hpsi_i(0);
		for (SizeType j = 0; j < dim; ++j)
			hpsi_i += denseH(i, j) * psi1[j];
		raw += std::conj(psi1[i]) * hpsi_i;
	}

	RealType    onsite    = 0;
	const auto& diagFixed = Acc::diagFixed(solver);
	for (SizeType m = 0; m < dim; ++m)
		onsite += diagFixed[m] * std::norm(psi1[m]);

	// Written out explicitly (not calling the production 0.5 factor) --
	// see this test's top comment for why.
	const RealType expectedEkin = RealType(0.5) * std::real(raw - ComplexType(onsite));
	const RealType gotEkin      = Acc::computeKineticEnergyGBEKSector(solver, 1);

	INFO("expected (dense) = " << expectedEkin << ", got (production CSR) = " << gotEkin);
	CHECK(gotEkin == Catch::Approx(expectedEkin).margin(1e-10));
	// Sanity: this Ekin must be nonzero, or the test would trivially pass
	// regardless of the 0.5 factor (as Test A above documents).
	CHECK(std::abs(gotEkin) > 1e-6);
}
