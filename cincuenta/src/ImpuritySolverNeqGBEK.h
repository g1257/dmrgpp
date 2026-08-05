#ifndef IMPURITYSOLVER_NEQ_GBEK_H
#define IMPURITYSOLVER_NEQ_GBEK_H

#include "CincuentaInputCheck.h"
#include "ImpuritySolverNeqBase.h"
#include "ImpuritySolverNeqExactDiag.h"
#include "KadanoffBaym.h"
#include "LanczosPlusPlus/src/Engine/DefaultSymmetry.h"
#include "LanczosPlusPlus/src/Engine/InputCheck.h"
#include "LanczosPlusPlus/src/Engine/InternalProductStored.h"
#include "LanczosPlusPlus/src/Engine/LabeledOperator.h"
#include "LanczosPlusPlus/src/Engine/LanczosGlobals.h"
#include "LanczosPlusPlus/src/Engine/ModelSelector.h"
#include "NeqBathDecomposition.h"
#include "ParamsNeqDmftSolver.h"
#include <PsimagLite/LanczosSolver.h>
#include <PsimagLite/Matrix.h>
#include <PsimagLite/ParametersForSolver.h>
#include <PsimagLite/PsimagLite.h>
#include <PsimagLite/Vector.h>

#include <PsimagLite/CrsMatrix.h>

#include <algorithm>
#include <cassert>
#include <fstream>
#include <iomanip>
#include <limits>
#include <memory>
#include <stdexcept>
#include <vector>

// Non-equilibrium impurity solver implementing the GBEK two-bath scheme:
//   Gramsch, Balzer, Eckstein, Kollar — PRB 88, 235106 (2013).
//
// The hybridization Δ is split as Δ = Δ⁻ + Δ⁺:
//   Δ⁻ : first bath — equilibrium memory from fixed {V_α, ε_α}
//   Δ⁺ : second bath — neq dynamics via rank-L Cholesky of i·Δ⁺_<
//
// For L=0 (first bath only), delegates entirely to ImpuritySolverNeqExactDiag
// (exact Lehmann representation), which is the correct L=0 limit.
//
// For L>0, the extended Fock space includes 2L second-bath sites:
//   L "empty" sites (initially unoccupied, ε=0) coupled to impurity via V^+_{n,p}
//   L "occupied" sites (initially doubly occupied, ε=0) coupled to impurity via V^+_{n,p}
// The many-body state is propagated step by step under the time-dependent H(t_n)
// using a Krylov Lanczos matrix exponential.
//
// Green's functions are computed as inner products (GBEK Eqs. 72-74):
//   G^<(t_j, t_n) = +i <Ψ(t_j) | Ψ(t_n)>   [N-1 sector]
//   G^>(t_n, t_j) = -i <Φ(t_n) | Φ(t_j)>   [N+1 sector]
// G^{Left}(t_n, τ) is reused from ImpuritySolverNeqExactDiag (Matsubara-based).
//
// Gα/Gβ spin-seed averaging (GBEK Eq. 70): whenever the base filling is
// spin-imbalanced (nup_ != ndown_ -- e.g. the NeqAtomicLimit single-atom
// seed nup=1, ndown=0), a single extended-Fock-space ground state is
// spin-polarized and does not represent the paramagnetic, particle-hole
// symmetric physics on its own. The paper's remedy: run the WHOLE
// extended-system calculation twice -- once seeded as if the impurity's
// extra electron were up (system alpha, nup_ext=nup+L, ndown_ext=ndown+L),
// once as if it were down (system beta, nup_ext=ndown+L, ndown_ext=nup+L)
// -- and average the resulting (always up-channel-probed) Green's
// functions. When nup_==ndown_ (e.g. the near-atomic path's half-filled
// multi-site bath), alpha and beta are the identical configuration, so
// this reduces to a no-op and only one calculation is actually run.
namespace Dmft {

template <typename ComplexOrRealType>
class ImpuritySolverNeqGBEK : public ImpuritySolverNeqBase<ComplexOrRealType> {

public:

	using BaseType          = ImpuritySolverNeqBase<ComplexOrRealType>;
	using RealType          = typename BaseType::RealType;
	using ComplexType       = typename BaseType::ComplexType;
	using VectorRealType    = typename BaseType::VectorRealType;
	using KBType            = typename BaseType::KBType;
	using InputNgType       = typename BaseType::InputNgType;
	using ParamsNeqType     = ParamsNeqDmftSolver<ComplexOrRealType>;
	using ExactDiagType     = ImpuritySolverNeqExactDiag<ComplexOrRealType>;
	using DecompType        = NeqBathDecomposition<ComplexOrRealType>;
	using VectorComplexType = typename PsimagLite::Vector<ComplexType>::Type;
	using MatrixComplexType = PsimagLite::Matrix<ComplexType>;
	using MatrixType        = PsimagLite::Matrix<ComplexOrRealType>;

	// LanczosPlusPlus / PsimagLite types (parallel to ImpuritySolverNeqExactDiag)
	using LppInputReadable =
	    typename PsimagLite::InputNg<LanczosPlusPlus::InputCheck>::Readable;
	using GeometryType = PsimagLite::
	    Geometry<ComplexOrRealType, LppInputReadable, LanczosPlusPlus::LanczosGlobals>;
	using ModelSelectorType
	    = LanczosPlusPlus::ModelSelector<ComplexOrRealType, GeometryType, LppInputReadable>;
	using ModelBaseType
	    = LanczosPlusPlus::LanczosModelBase<ComplexOrRealType, GeometryType, LppInputReadable>;
	using BasisBaseType       = typename ModelBaseType::BasisBaseType;
	using DefaultSymmetryType = LanczosPlusPlus::DefaultSymmetry<GeometryType, BasisBaseType>;
	using InternalProductStoredType
	    = LanczosPlusPlus::InternalProductStored<ModelBaseType, DefaultSymmetryType>;
	using LabeledOperatorType    = LanczosPlusPlus::LabeledOperator;
	using WordType               = LanczosPlusPlus::LanczosGlobals::WordType;
	using PairIntType            = LanczosPlusPlus::LanczosGlobals::PairIntType;
	using LanczosSolverForGSType = PsimagLite::LanczosSolver<InternalProductStoredType>;
	using CrsMatrixComplexType   = PsimagLite::CrsMatrix<ComplexType>;

	ImpuritySolverNeqGBEK(const ParamsNeqType& params, typename InputNgType::Readable& io)
	    : bathRank_(params.neqBathRank)
	    , params_(params)
	    , exactDiag_(params, io)
	    , nup_(0)
	    , ndown_(0)
	    , nBath_(0)
	    , nsites_ext_(0)
	{
		if (bathRank_ > 0) {
			io.readline(nup_, "TargetElectronsUp=");
			io.readline(ndown_, "TargetElectronsDown=");
		}
	}

	// Initialise the impurity solver and bath decomposition.
	// For L=0: delegates to ExactDiag (Lehmann representation).
	// For L>0: also builds the extended Fock space and seeds time-propagated states.
	void solve(const VectorRealType& bathParams) override
	{
		const RealType beta = params_.eqParams.ficticiousBeta;
		const RealType mu   = 0;

		decomp_ = std::make_unique<DecompType>(
		    bathRank_,
		    beta,
		    mu,
		    bathParams,
		    params_.nT,
		    params_.eqParams.nMatsubaras,
		    params_.dt,
		    params_.eqParams.ficticiousBeta
		        / static_cast<RealType>(params_.eqParams.nMatsubaras));

		exactDiag_.solve(bathParams);

		if (bathRank_ > 0)
			solveLplus(bathParams);
	}

	// Advance the Cholesky decomposition to step n and (for L>0) invalidate
	// the propagated state so the next computeGimp call re-propagates.
	void prepareTimeStep(int n, const KBType& delta) override
	{
		if (decomp_)
			decomp_->update(n, delta);
		if (bathRank_ > 0 && n > 0) {
			// V[n] just changed; invalidate cached PsiHist[n] and beyond, for
			// both spin configurations.
			sectorAlpha_.propagatedThrough
			    = std::min(sectorAlpha_.propagatedThrough, n - 1);
			if (!sameConfig_)
				sectorBeta_.propagatedThrough
				    = std::min(sectorBeta_.propagatedThrough, n - 1);
		}
	}

	// Fill G_imp KB components for all (t_n, t_j) with j <= n.
	// For L=0: delegates to ExactDiag (Lehmann).
	// For L>0: ExactDiag fills G^{Left}; GBEK inner products overwrite G^< and G^R.
	void computeGimp(KBType& gimp, int n) const override
	{
		// Always call ExactDiag for G^{Left} and Matsubara components
		exactDiag_.computeGimp(gimp, n);
		if (bathRank_ == 0)
			return;

		// Ensure states are propagated up to step n (lazy, idempotent)
		ensurePropagated(n);

		// Double occupation / kinetic energy at this step (paper Figs. 9-10)
		// -- cheap next to the propagation/Green's-function work above (one
		// extra sparse matvec plus an O(dim) diagonal loop; see
		// computeKineticEnergyGBEKSector).
		doccHistory_[static_cast<SizeType>(n)] = computeDoccGBEK(n);
		ekinHistory_[static_cast<SizeType>(n)] = computeKineticEnergyGBEK(n);

		// Overwrite G^< and G^R with GBEK inner-product formulas. Both rows
		// are computed in one backward Krylov sweep each (see
		// gLesserRowGBEKSector/gGreaterRowGBEKSector), not per-(n,j).
		const std::vector<ComplexType> gLessRow = gLesserRowGBEK(n);
		const std::vector<ComplexType> gGrtrRow = gGreaterRowGBEK(n);
		for (int j = 0; j <= n; ++j) {
			const SizeType    sn    = static_cast<SizeType>(n);
			const SizeType    sj    = static_cast<SizeType>(j);
			const ComplexType gLess = gLessRow[static_cast<SizeType>(j)];
			const ComplexType gGrtr = gGrtrRow[static_cast<SizeType>(j)];
			gimp.lesser(sn, sj)     = gLess;
			gimp.retarded(sn, sj)   = gGrtr - gLess;
		}
		// Anti-Hermitian: G^<(t_j, t_n) = -conj(G^<(t_n, t_j)) for j < n
		for (int j = 0; j < n; ++j)
			gimp.lesser(static_cast<SizeType>(j), static_cast<SizeType>(n))
			    = -std::conj(
			        gimp.lesser(static_cast<SizeType>(n), static_cast<SizeType>(j)));
	}

	// Matsubara components come from ExactDiag in both L=0 and L>0 cases.
	const KBType& gimp() const override { return exactDiag_.gimp(); }

	const DecompType* decomposition() const { return decomp_.get(); }

	void dumpPlusBath(const std::string& filename) const override
	{
		if (decomp_)
			decomp_->dumpPlusBath(filename);
	}

	// Write one line per time step actually visited (n=0..propagatedThrough):
	// "t docc Ekin Eint Etot", fixed precision, mirroring KadanoffBaym::dump's
	// plain-text style. Eint/Etot are trivial derived quantities
	// (Eint=U*(docc-1/4), Etot=Ekin+Eint); docc/Ekin come from doccHistory_/
	// ekinHistory_, filled progressively by computeGimp. See paper Figs. 9-10
	// and project_gbek_ekin_factor2_bugfix for the 0.5 factor in Ekin.
	void dumpDoccAndEnergy(const std::string& filename) const override
	{
		if (bathRank_ == 0)
			return;
		std::ofstream fout(filename.c_str());
		const int     nMax = sectorAlpha_.propagatedThrough;
		for (int n = 0; n <= nMax; ++n) {
			const SizeType sn   = static_cast<SizeType>(n);
			const RealType t    = static_cast<RealType>(n) * params_.dt;
			const RealType docc = doccHistory_[sn];
			const RealType ekin = ekinHistory_[sn];
			const RealType eint = params_.uFinal * (docc - RealType(0.25));
			const RealType etot = ekin + eint;
			fout << std::fixed << std::setprecision(10) << t << " " << docc << " "
			     << ekin << " " << eint << " " << etot << "\n";
		}
	}

	// Raw Cholesky factor V_(n,p), for row-by-row comparison against an
	// independent offline trace (see NeqBathDecomposition::dumpV).
	void dumpV(const std::string& filename) const override
	{
		if (decomp_)
			decomp_->dumpV(filename);
	}

private:

	// Position and update info for one second-bath variable entry in the CSR matrix.
	struct VarEntry {
		int      nnzIdx; // index into CrsMatrix values array (matches setValues arg type)
		SizeType p; // second-bath rank index
		bool     isConj; // true → store conj(vMid[p])*sign, false → vMid[p]*sign
		int      sign; // Jordan-Wigner sign (±1)
	};

	// ========== Gα/Gβ per-configuration extended-Fock state ==========
	//
	// Everything that depends on WHICH spin the impurity's extra electron
	// (nup_ext,ndown_ext) is seeded with. System alpha uses
	// (nup+L, ndown+L) (impurity's extra electron is up); system beta
	// swaps to (ndown+L, nup+L) (as if it were down). Both are always
	// probed with the up-spin operator (see seedState/gLesserGBEKSector),
	// matching GBEK Eq. 70's Gα_σ/Gβ_σ for σ=up.
	// PhiNHist: the post-quench N-sector reference trajectory (forward-
	// propagated pre-quench GS). bStates[n]/dStates[n]: c_{imp,up}/c†_{imp,up}
	// applied to PhiNHist[n] AT EACH n -- this reseeding is required for
	// correctness (see propagateOneStep's doc comment); it replaces an
	// earlier, incorrect scheme that seeded once at n=0 and propagated the
	// (N∓1)-sector state forward using only the (N∓1)-sector Hamiltonian.
	struct ExtendedSector {
		SizeType                               dim1Nm1 = 0, dim1Np1 = 0, dim1N = 0;
		std::vector<WordType>                  upWordsNm1, dnWordsNm1;
		std::vector<WordType>                  upWordsNp1, dnWordsNp1;
		std::vector<WordType>                  upWordsN, dnWordsN;
		CrsMatrixComplexType                   csrNm1, csrNp1, csrN;
		std::vector<VarEntry>                  varNm1, varNp1, varN;
		CrsMatrixComplexType                   cUpNm1, cUpDagNp1;
		mutable std::vector<VectorComplexType> PhiNHist, bStates, dStates;
		mutable int                            propagatedThrough = -1;
		// Fixed onsite (U + potential) diagonal of the N-sector Hamiltonian,
		// indexed by N-sector basis index -- captured once in buildHextCSR
		// (it never changes) so <Ekin(t)> can be read off as
		// <PhiNHist[n]|csrN|PhiNHist[n]> minus this term, without needing to
		// re-derive or extract the diagonal from the CSR matrix later. See
		// computeKineticEnergyGBEKSector.
		std::vector<RealType> diagFixed;
	};

	// ========== L>0 setup ==========

	// Build extended Fock space, diagonalise pre-quench H, seed PsiHist/PhiHist
	// for both spin configurations (alpha always; beta only if nup_ != ndown_,
	// since otherwise it is identical to alpha -- see class docstring).
	void solveLplus(const VectorRealType& bathParams)
	{
		const SizeType L = bathRank_;
		nBath_           = bathParams.size() / 2;
		nsites_ext_      = nBath_ + 1 + 2 * L;

		firstBathHop_.resize(nBath_);
		firstBathEps_.resize(nBath_);
		for (SizeType i = 0; i < nBath_; ++i) {
			firstBathHop_[i] = bathParams[i];
			firstBathEps_[i] = bathParams[nBath_ + i];
		}

		// Post-quench on-site potentials (ε=0 for second bath)
		potPost_.assign(nsites_ext_, RealType(0));
		potPost_[0] = -RealType(0.5) * params_.uFinal;
		for (SizeType i = 0; i < nBath_; ++i)
			potPost_[i + 1] = firstBathEps_[i];

		sameConfig_ = (nup_ == ndown_);

		doccHistory_.assign(params_.nT + 1, RealType(0));
		ekinHistory_.assign(params_.nT + 1, RealType(0));

		buildSector(sectorAlpha_, nup_ + L, ndown_ + L, bathParams);
		if (!sameConfig_)
			buildSector(sectorBeta_, ndown_ + L, nup_ + L, bathParams);
	}

	// Build one spin configuration's extended Fock space, diagonalise its
	// pre-quench Hamiltonian, and seed PsiHist[0]/PhiHist[0] -- shared logic
	// for both system alpha and system beta (GBEK Eq. 70).
	void buildSector(ExtendedSector&       sector,
	                 SizeType              nupExt,
	                 SizeType              ndownExt,
	                 const VectorRealType& bathParams)
	{
		const SizeType L = bathRank_;

		// Extended hoppings (real; second bath starts at 0, updated later via Cholesky)
		VectorRealType hopExt(nBath_ + 2 * L, RealType(0));
		for (SizeType i = 0; i < nBath_; ++i)
			hopExt[i] = firstBathHop_[i];

		// Pre-quench extended potential: large ±ε for second bath to fix initial occupation
		const RealType bigEps = RealType(500);
		VectorRealType potPre(nsites_ext_, RealType(0));
		potPre[0] = -RealType(0.5) * params_.uInitial;
		for (SizeType i = 0; i < nBath_; ++i)
			potPre[i + 1] = firstBathEps_[i];
		for (SizeType p = 0; p < L; ++p)
			potPre[nBath_ + 1 + p] = +bigEps; // empty
		for (SizeType p = 0; p < L; ++p)
			potPre[nBath_ + 1 + L + p] = -bigEps; // occupied

		// Build extended pre-quench model
		const std::string inputPre = buildLanczosInput(
		    params_.uInitial, nupExt, ndownExt, hopExt, potPre, nsites_ext_);

		LanczosPlusPlus::InputCheck                                          ic;
		typename PsimagLite::InputNg<LanczosPlusPlus::InputCheck>::Writeable ioW(ic,
		                                                                         inputPre);
		LppInputReadable                                                     ioR(ioW);
		GeometryType                                                         geom(ioR);
		ModelSelectorType                                                    ms(ioR, geom);
		const ModelBaseType&                                                 model = ms();

		// Diagonalise N sector → pre-quench GS.
		// Full diag is only feasible up to nsites_ext_<=8 (dense, O(dim^3)/O(dim^2)
		// time/memory -- infeasible well before L=5, dim~213k); for larger
		// systems use Lanczos GS, seeded with the KNOWN atomic-limit product
		// state (see lanczosGS's doc comment) rather than a naive uniform
		// vector -- see project memory project_gbek_cpp_lanczosGS_bug for the
		// bug this fixes (uniform-vector Lanczos converged to the WRONG
		// extremal state at L=4: reported energy -2807.96 vs. the true -4000).
		VectorRealType energiesN;
		MatrixType     eigvecsN;
		if (nsites_ext_ <= 8) {
			diagWithBasis(model, model.basis(), geom, energiesN, eigvecsN);
		} else {
			// Impurity's own occupation within this sector: each of the L
			// "occ" second-bath sites contributes exactly 1 up + 1 down
			// electron (see potPre's bigEps assignments just above), so
			// whatever nupExt/ndownExt has left over is what the impurity
			// itself must carry (0 or 1 each, enforced by the atomic-limit
			// nup_+ndown_==1 constraint upstream).
			assert(nupExt >= L && ndownExt >= L);
			const SizeType impurityUp   = nupExt - L;
			const SizeType impurityDown = ndownExt - L;
			eigvecsN = lanczosGS(model, geom, nBath_, L, impurityUp, impurityDown);
		}

		// Build N-1 and N+1 sector bases for the extended system.
		// model owns the returned basis pointers via its internal garbage_
		// list; do not wrap in unique_ptr (double-delete crashes
		// ~HubbardOneOrbital) -- same footgun documented in
		// ImpuritySolverNeqExactDiag.h's solve().
		const BasisBaseType& bNm1 = *model.createBasis(nupExt - 1, ndownExt);
		const BasisBaseType& bNp1 = *model.createBasis(nupExt + 1, ndownExt);

		// Build fast-lookup tables from sorted Fock words
		buildFockLookup(bNm1, sector.upWordsNm1, sector.dnWordsNm1, sector.dim1Nm1);
		buildFockLookup(bNp1, sector.upWordsNp1, sector.dnWordsNp1, sector.dim1Np1);
		buildFockLookup(model.basis(), sector.upWordsN, sector.dnWordsN, sector.dim1N);

		// Build CSR sparse matrices for H_ext and validate against applyHext
		buildHextCSR(sector.upWordsNm1,
		             sector.dnWordsNm1,
		             sector.dim1Nm1,
		             sector.csrNm1,
		             sector.varNm1);
		buildHextCSR(sector.upWordsNp1,
		             sector.dnWordsNp1,
		             sector.dim1Np1,
		             sector.csrNp1,
		             sector.varNp1);
		// N-sector Hamiltonian: needed to propagate the reference trajectory
		// PhiNHist forward (see propagateOneStep's doc comment). Also capture
		// its fixed onsite diagonal into sector.diagFixed, for
		// computeKineticEnergyGBEKSector.
		buildHextCSR(sector.upWordsN,
		             sector.dnWordsN,
		             sector.dim1N,
		             sector.csrN,
		             sector.varN,
		             &sector.diagFixed);
		checkApplyHextEquivalence(sector);

		// Sparse c_{imp,↑} / c†_{imp,↑} operators, N-sector -> N-1/N+1 sector.
		// Built once here (basisN/basisNm1/basisNp1 still in scope) and reused
		// at every time step to reseed against PhiNHist[n] -- see
		// propagateOneStep's doc comment for why this reseeding matters.
		sector.cUpNm1 = buildCOperatorCSR(
		    model.basis(), bNm1, LabeledOperatorType::Label::OPERATOR_C);
		sector.cUpDagNp1 = buildCOperatorCSR(
		    model.basis(), bNp1, LabeledOperatorType::Label::OPERATOR_CDAGGER);

		// Allocate history arrays
		const SizeType nSteps = params_.nT + 1;
		const SizeType dimN   = model.basis().size();
		sector.PhiNHist.assign(nSteps, VectorComplexType(dimN, ComplexType(0)));
		sector.bStates.assign(nSteps, VectorComplexType(bNm1.size(), ComplexType(0)));
		sector.dStates.assign(nSteps, VectorComplexType(bNp1.size(), ComplexType(0)));

		// Seed step 0: PhiNHist[0] = pre-quench GS; bStates[0]/dStates[0]
		// follow by applying the (now reusable) c/c† operators to it -- must
		// reproduce the old one-off PsiHist[0]/PhiHist[0] construction exactly.
		for (SizeType m = 0; m < dimN; ++m)
			sector.PhiNHist[0][m] = ComplexType(eigvecsN(m, 0));
		sparseMatVec(sector.cUpNm1, sector.PhiNHist[0], sector.bStates[0]);
		sparseMatVec(sector.cUpDagNp1, sector.PhiNHist[0], sector.dStates[0]);
		sector.propagatedThrough = 0;
	}

	// Build the sparse CSR for c_{imp,↑} (opLabel=OPERATOR_C) or c†_{imp,↑}
	// (opLabel=OPERATOR_CDAGGER), mapping the N-sector Fock basis (dimN
	// columns) to a companion sector (N-1 or N+1, dimTo rows). Always the
	// up-spin operator, regardless of whether this sector is system alpha or
	// system beta -- see class docstring. Applicable to ANY N-sector state
	// vector, not just the pre-quench GS (unlike the old seedState, which
	// only ever built the t=0 seed).
	static CrsMatrixComplexType buildCOperatorCSR(const BasisBaseType&                basisN,
	                                              const BasisBaseType&                basisTo,
	                                              typename LabeledOperatorType::Label opLabel)
	{
		const SizeType            dimN    = basisN.size();
		const SizeType            dimTo   = basisTo.size();
		const SizeType            spinUp  = LanczosPlusPlus::LanczosGlobals::SPIN_UP;
		const SizeType            spinDn  = LanczosPlusPlus::LanczosGlobals::SPIN_DOWN;
		const SizeType            impSite = 0;
		const LabeledOperatorType op(opLabel);

		struct Entry {
			SizeType    row, col;
			ComplexType val;
		};
		std::vector<Entry> entries;
		entries.reserve(dimN);

		for (SizeType m = 0; m < dimN; ++m) {
			const WordType    ket1 = basisN(m, spinUp);
			const WordType    ket2 = basisN(m, spinDn);
			const PairIntType bra
			    = basisTo.getBraIndex(ket1, ket2, op, impSite, spinUp, 0);
			if (bra.first < 0)
				continue;
			const int sign = basisN.doSignGf(ket1, ket2, impSite, spinUp, 0);
			entries.push_back(
			    { static_cast<SizeType>(bra.first), m, ComplexType(sign) });
		}

		std::sort(entries.begin(),
		          entries.end(),
		          [](const Entry& a, const Entry& b)
		          { return a.row < b.row || (a.row == b.row && a.col < b.col); });

		CrsMatrixComplexType csr;
		csr.resize(dimTo, dimN);
		csr.reserve(entries.size());
		SizeType idx = 0;
		int      nnz = 0;
		for (SizeType row = 0; row < dimTo; ++row) {
			csr.setRow(row, nnz);
			while (idx < entries.size() && entries[idx].row == row) {
				csr.pushCol(entries[idx].col);
				csr.pushValue(entries[idx].val);
				++nnz;
				++idx;
			}
		}
		csr.setRow(dimTo, nnz);
		return csr;
	}

	// Extract sorted up/dn Fock words from a basis for O(log n) index lookup.
	static void buildFockLookup(const BasisBaseType&   basis,
	                            std::vector<WordType>& upWords,
	                            std::vector<WordType>& dnWords,
	                            SizeType&              dim1)
	{
		const SizeType dim    = basis.size();
		const SizeType spinUp = LanczosPlusPlus::LanczosGlobals::SPIN_UP;
		const SizeType spinDn = LanczosPlusPlus::LanczosGlobals::SPIN_DOWN;

		if (dim == 0) {
			dim1 = 0;
			return;
		}

		// Determine dim1: period of the up-spin word (it cycles with period = #distinct up
		// words)
		const WordType firstUp = basis(0, spinUp);
		dim1                   = 1;
		while (dim1 < dim && basis(dim1, spinUp) != firstUp)
			++dim1;

		upWords.resize(dim1);
		for (SizeType x = 0; x < dim1; ++x)
			upWords[x] = basis(x, spinUp);

		const SizeType dim2 = dim / dim1;
		dnWords.resize(dim2);
		for (SizeType y = 0; y < dim2; ++y)
			dnWords[y] = basis(y * dim1, spinDn);
	}

	// ========== Precomputed CSR sparse Hamiltonian ==========

	// Build the CSR sparse matrix for H_ext in one sector.
	// Fixed entries (diagonal + first-bath hops) go directly into matrix values.
	// Variable entries (second-bath hops, time-dependent) are stored as zero initially;
	// their positions are recorded in varEntries for O(1) update at each time step.
	// Asserts that no (row,col) slot is simultaneously fixed and variable.
	void buildHextCSR(const std::vector<WordType>& upWords,
	                  const std::vector<WordType>& dnWords,
	                  SizeType                     dim1,
	                  CrsMatrixComplexType&        csr,
	                  std::vector<VarEntry>&       varEntries,
	                  std::vector<RealType>*       diagOut = nullptr) const
	{
		const SizeType dim2   = dnWords.size();
		const SizeType dim    = dim1 * dim2;
		const SizeType spinUp = LanczosPlusPlus::LanczosGlobals::SPIN_UP;
		if (diagOut)
			diagOut->assign(dim, RealType(0));

		struct TripletEntry {
			SizeType    row, col;
			ComplexType fixedVal;
			bool        isVar;
			SizeType    p;
			bool        isConj;
			int         sign;
		};

		std::vector<TripletEntry> triplets;
		triplets.reserve(dim * (1 + 2 * (2 * nBath_ + 4 * bathRank_)));

		for (SizeType m = 0; m < dim; ++m) {
			const SizeType y    = m / dim1;
			const SizeType x    = m % dim1;
			const WordType up_m = upWords[x];
			const WordType dn_m = dnWords[y];

			// Diagonal term
			RealType diag = 0;
			if ((up_m & WordType(1)) && (dn_m & WordType(1)))
				diag += params_.uFinal;
			for (SizeType i = 0; i < nsites_ext_; ++i) {
				if ((up_m >> i) & WordType(1))
					diag += potPost_[i];
				if ((dn_m >> i) & WordType(1))
					diag += potPost_[i];
			}
			if (diag != RealType(0))
				triplets.push_back({ m, m, ComplexType(diag), false, 0, false, 1 });
			if (diagOut)
				(*diagOut)[m] = diag;

			for (SizeType spin = 0; spin < 2; ++spin) {
				const WordType w       = (spin == spinUp) ? up_m : dn_m;
				const WordType w_other = (spin == spinUp) ? dn_m : up_m;

				auto tryHop = [&](SizeType    from,
				                  SizeType    to,
				                  bool        isVar,
				                  SizeType    p,
				                  bool        isConj,
				                  ComplexType fixedT)
				{
					if (!((w >> from) & WordType(1)))
						return;
					if ((w >> to) & WordType(1))
						return;
					const WordType new_w
					    = (w & ~(WordType(1) << from)) | (WordType(1) << to);
					const WordType new_up = (spin == spinUp) ? new_w : w_other;
					const WordType new_dn = (spin == spinUp) ? w_other : new_w;
					const SizeType m_new  = lookupIndex(
                                            new_up, new_dn, upWords, dnWords, dim1, dim2);
					if (m_new == std::numeric_limits<SizeType>::max())
						return;
					const int         sign = computeJWSign(w, from, to);
					const ComplexType val
					    = isVar ? ComplexType(0) : fixedT * ComplexType(sign);
					triplets.push_back(
					    { m_new, m, val, isVar, p, isConj, sign });
				};

				// First bath (fixed, real hoppings V_α)
				for (SizeType a = 0; a < nBath_; ++a) {
					const ComplexType Va(firstBathHop_[a]);
					tryHop(0, a + 1, false, 0, false, Va);
					tryHop(a + 1, 0, false, 0, false, Va);
				}

				// Second bath (variable: hop ↔ empty-site, ↔ occupied-site)
				for (SizeType p = 0; p < bathRank_; ++p) {
					const SizeType emptyS = 1 + nBath_ + p;
					const SizeType occS   = 1 + nBath_ + bathRank_ + p;
					tryHop(0, emptyS, true, p, false, ComplexType(0));
					tryHop(emptyS, 0, true, p, true, ComplexType(0));
					// Occupied-orbital coupling has the OPPOSITE conjugation
					// from the empty orbital's (GBEK Eq. 49/52a: -iLambda^<_+
					// couples the occupied orbital via V(t)V(t')^*, while
					// +iLambda^>_+ =
					// (-iLambda^<_+)^* couples the empty orbital -- see
					// gbek_selfconsistency.py::build_templates's docstring,
					// which documents this exact asymmetry and was validated
					// against the paper directly). c^dag_occ c_imp must be
					// conj(V), and c^dag_imp c_occ must be V -- opposite of
					// what was here before.
					tryHop(0, occS, true, p, true, ComplexType(0));
					tryHop(occS, 0, true, p, false, ComplexType(0));
				}
			}
		}

		// Sort by (row, col) so we can build row-major CSR.
		std::sort(triplets.begin(),
		          triplets.end(),
		          [](const TripletEntry& a, const TripletEntry& b)
		          { return a.row < b.row || (a.row == b.row && a.col < b.col); });

		// Assert no fixed/var collision at the same (row,col) slot.
		for (SizeType i = 1; i < triplets.size(); ++i) {
			if (triplets[i].row == triplets[i - 1].row
			    && triplets[i].col == triplets[i - 1].col)
				assert(false
				       && "buildHextCSR: fixed and variable hop collide at same "
				          "(row,col)");
		}

		varEntries.clear();
		csr.resize(dim, dim);
		csr.reserve(triplets.size());

		int      nnzCounter = 0;
		SizeType tripletIdx = 0;
		for (SizeType row = 0; row < dim; ++row) {
			csr.setRow(row, nnzCounter);
			while (tripletIdx < triplets.size() && triplets[tripletIdx].row == row) {
				const TripletEntry& t = triplets[tripletIdx];
				csr.pushCol(t.col);
				csr.pushValue(t.fixedVal);
				if (t.isVar)
					varEntries.push_back({ nnzCounter, t.p, t.isConj, t.sign });
				++nnzCounter;
				++tripletIdx;
			}
		}
		csr.setRow(dim, nnzCounter);
	}

	// Overwrite second-bath (variable) CSR entries with current midpoint hoppings.
	// All other entries retain their fixed values from buildHextCSR.
	void updateCSR(CrsMatrixComplexType&           csr,
	               const std::vector<VarEntry>&    varEntries,
	               const std::vector<ComplexType>& vMid) const
	{
		for (const VarEntry& ve : varEntries) {
			const ComplexType Vp = ve.isConj ? std::conj(vMid[ve.p]) : vMid[ve.p];
			csr.setValues(ve.nnzIdx, Vp * ComplexType(ve.sign));
		}
	}

	// ========== Time propagation ==========

	// Ensure PsiHist_[0..n] and PhiHist_[0..n] are populated with the current V values.
	void ensurePropagated(int n) const
	{
		ensureSectorPropagated(sectorAlpha_, n);
		if (!sameConfig_)
			ensureSectorPropagated(sectorBeta_, n);
	}

	void ensureSectorPropagated(ExtendedSector& sector, int n) const
	{
		for (int k = sector.propagatedThrough + 1; k <= n; ++k)
			propagateOneStep(sector, k);
		if (n > sector.propagatedThrough)
			sector.propagatedThrough = n;
	}

	// Midpoint second-bath hoppings for the interval ending at step n:
	// V_mid = (V[n-1] + V[n]) / 2. Shared by forward propagation of PhiNHist
	// and by the backward sweeps in gLesserRowGBEKSector/gGreaterRowGBEKSector
	// -- the same physical Hamiltonian governs interval [t_{n-1}, t_n]
	// regardless of which direction it is traversed.
	std::vector<ComplexType> computeVMid(int n) const
	{
		std::vector<ComplexType> vMid(bathRank_);
		for (SizeType p = 0; p < bathRank_; ++p) {
			const ComplexType vPrev = decomp_->Vplus(n - 1, static_cast<int>(p));
			const ComplexType vCurr = decomp_->Vplus(n, static_cast<int>(p));
			vMid[p]                 = RealType(0.5) * (vPrev + vCurr);
		}
		return vMid;
	}

	// Advance the N-sector reference trajectory PhiNHist from step n-1 to n
	// under the midpoint Hamiltonian, then reseed c_{imp,up}/c†_{imp,up}
	// against PhiNHist[n] to get bStates[n]/dStates[n].
	//
	// This reseeding-at-every-step is required for correctness: the two-time
	// Green's function is G^<(t,t') = i<c†(t')c(t)>, and in the Heisenberg
	// picture c(t)|ψ0> = U_{N-1}(0,t) · c · U_N(t,0)|ψ0> -- i.e. c must be
	// applied to the N-sector state AT TIME t, not at t=0. An earlier version
	// of this solver instead applied c ONCE (at t=0) and propagated the
	// resulting (N-1)-sector state forward using only the (N-1)-sector
	// Hamiltonian; that shortcut is valid only if c commutes with H, which
	// fails here because of the Hubbard U term ([U n_up n_dn, c_up] != 0).
	// The bug showed up as a systematically suppressed off-diagonal
	// (two-time) imaginary part, growing with |n-j| while sparing the
	// diagonal -- confirmed against an independent brute-force
	// (dense-matrix-exponential) reconstruction in
	// cincuenta/TestSuite/gbek_reference/cross_check_seed_scheme.py.
	void propagateOneStep(ExtendedSector& sector, int n) const
	{
		assert(n > 0 && n < static_cast<int>(sector.PhiNHist.size()));

		const std::vector<ComplexType> vMid = computeVMid(n);

		updateCSR(sector.csrN, sector.varN, vMid);
		sector.PhiNHist[n]
		    = krylovExpmvCSR(sector.PhiNHist[n - 1], sector.csrN, params_.dt);

		sparseMatVec(sector.cUpNm1, sector.PhiNHist[n], sector.bStates[n]);
		sparseMatVec(sector.cUpDagNp1, sector.PhiNHist[n], sector.dStates[n]);
	}

	// Rectangular-safe sparse mat-vec: x = A*y. PsimagLite::CrsMatrix's own
	// matrixVectorProduct assumes a square matrix (it asserts
	// x.size()==y.size()); cUpNm1_/cUpDagNp1_ map between differently-sized
	// sectors, so they need this instead.
	static void sparseMatVec(const CrsMatrixComplexType& A,
	                         const VectorComplexType&    y,
	                         VectorComplexType&          x)
	{
		x.assign(A.rows(), ComplexType(0));
		for (SizeType i = 0; i < A.rows(); ++i)
			for (int k = A.getRowPtr(i); k < A.getRowPtr(i + 1); ++k)
				x[i] += A.getValue(static_cast<SizeType>(k))
				    * y[static_cast<SizeType>(A.getCol(static_cast<SizeType>(k)))];
	}

	// Krylov (Lanczos) approximation to exp(-i H dt) |psi⟩.
	// Uses the Lanczos recursion to project H into a small Krylov subspace,
	// then computes the matrix exponential there.
	VectorComplexType krylovExpmv(const VectorComplexType&        psi,
	                              const std::vector<ComplexType>& vMid,
	                              const std::vector<WordType>&    upWords,
	                              const std::vector<WordType>&    dnWords,
	                              SizeType                        dim1,
	                              RealType                        dt,
	                              int                             m = 40) const
	{
		const SizeType dim   = psi.size();
		const RealType norm0 = std::sqrt(std::real(innerProduct(psi, psi)));
		if (norm0 < RealType(1e-14))
			return psi;

		// Lanczos iteration to build K_m(H, psi/norm0)
		std::vector<VectorComplexType> Q;
		Q.reserve(m + 1);
		Q.push_back(psi);
		scaleVec(Q[0], ComplexType(RealType(1) / norm0));

		std::vector<RealType> alpha, beta;
		VectorComplexType     Hq(dim);
		int                   mUsed = 0;

		for (int j = 0; j < m; ++j) {
			applyHext(Q[j], Hq, vMid, upWords, dnWords, dim1);

			alpha.push_back(std::real(innerProduct(Q[j], Hq)));
			VectorComplexType w = Hq;
			axpy(w, -ComplexType(alpha[j]), Q[j]);
			if (j > 0)
				axpy(w, -ComplexType(beta[j - 1]), Q[j - 1]);

			const RealType b = std::sqrt(std::real(innerProduct(w, w)));
			mUsed            = j + 1;
			if (b < RealType(1e-10)) {
				std::cout << "  krylov early-exit at j=" << j << " / " << m << "\n";
				break; // invariant subspace
			}

			beta.push_back(b);
			Q.push_back(w);
			scaleVec(Q.back(), ComplexType(RealType(1) / b));
		}

		// Build complex tridiagonal T (real symmetric, stored as complex for
		// PsimagLite::diag)
		MatrixComplexType Tc(mUsed, mUsed, ComplexType(0));
		for (int j = 0; j < mUsed; ++j)
			Tc(j, j) = ComplexType(alpha[j]);
		for (int j = 0; j < mUsed - 1; ++j) {
			Tc(j + 1, j) = ComplexType(beta[j]);
			Tc(j, j + 1) = ComplexType(beta[j]);
		}

		// Diagonalise T_m: Tc gets overwritten with eigenvectors (columns)
		VectorRealType eigsTm;
		PsimagLite::diag(Tc, eigsTm, 'V');

		// Compute exp(-i T_m dt) e_1 in the eigenbasis and transform back
		// expCoeff[j] = Σ_k conj(Tc[0,k]) * exp(-i λ_k dt) * Tc[j,k]
		VectorComplexType expCoeff(mUsed, ComplexType(0));
		for (int k = 0; k < mUsed; ++k) {
			const ComplexType phase  = std::exp(ComplexType(0, -eigsTm[k] * dt));
			const ComplexType factor = std::conj(Tc(0, k)) * phase;
			for (int j = 0; j < mUsed; ++j)
				expCoeff[j] += factor * Tc(j, k);
		}

		// Result: norm0 * Σ_j expCoeff[j] * Q[j]
		VectorComplexType result(dim, ComplexType(0));
		for (int j = 0; j < mUsed; ++j) {
			const ComplexType c = norm0 * expCoeff[j];
			for (SizeType i = 0; i < dim; ++i)
				result[i] += c * Q[j][i];
		}
		return result;
	}

	// Krylov approximation using precomputed CSR sparse matrix.
	// Caller must call updateCSR before this to load current vMid values.
	// Identical algorithm to krylovExpmv; SpMV replaces applyHext.
	//
	// For dim small enough that the Krylov recursion below (unreorthogonalized
	// Lanczos, standard for this kind of propagator) nearly exhausts the
	// invariant subspace before its early-exit check fires, the recursion
	// becomes ill-conditioned in the classic Lanczos "ghost eigenvalue" sense:
	// ULP-level input differences (BLAS threading order, SIMD codepath) can
	// flip whether one more, barely-linearly-independent Krylov vector gets
	// built, changing the propagated result by ~1e-4 -- found via this class's
	// own seed-scheme Catch2 test becoming flaky (nBath=0 sectors, dim=3-9).
	// Below this threshold, go via exact dense diagonalization instead: cost
	// is O(dim^3), negligible at this size, and exact regardless of any
	// degenerate/near-degenerate eigenvalues.
	VectorComplexType krylovExpmvCSR(const VectorComplexType&    psi,
	                                 const CrsMatrixComplexType& csr,
	                                 RealType                    dt,
	                                 int                         m = 40) const
	{
		const SizeType dim   = psi.size();
		const RealType norm0 = std::sqrt(std::real(innerProduct(psi, psi)));
		if (norm0 < RealType(1e-14))
			return psi;

		if (dim <= 64)
			return denseExpmv(csr, psi, dt);

		std::vector<VectorComplexType> Q;
		Q.reserve(m + 1);
		Q.push_back(psi);
		scaleVec(Q[0], ComplexType(RealType(1) / norm0));

		std::vector<RealType> alpha, beta;
		VectorComplexType     Hq(dim, ComplexType(0));
		int                   mUsed = 0;

		for (int j = 0; j < m; ++j) {
			std::fill(Hq.begin(), Hq.end(), ComplexType(0));
			csr.matrixVectorProduct(Hq, Q[j]);

			alpha.push_back(std::real(innerProduct(Q[j], Hq)));
			VectorComplexType w = Hq;
			axpy(w, -ComplexType(alpha[j]), Q[j]);
			if (j > 0)
				axpy(w, -ComplexType(beta[j - 1]), Q[j - 1]);

			const RealType b = std::sqrt(std::real(innerProduct(w, w)));
			mUsed            = j + 1;
			if (b < RealType(1e-10)) {
				std::cout << "  krylov CSR early-exit at j=" << j << " / " << m
				          << "\n";
				break;
			}

			beta.push_back(b);
			Q.push_back(w);
			scaleVec(Q.back(), ComplexType(RealType(1) / b));
		}

		MatrixComplexType Tc(mUsed, mUsed, ComplexType(0));
		for (int j = 0; j < mUsed; ++j)
			Tc(j, j) = ComplexType(alpha[j]);
		for (int j = 0; j < mUsed - 1; ++j) {
			Tc(j + 1, j) = ComplexType(beta[j]);
			Tc(j, j + 1) = ComplexType(beta[j]);
		}

		VectorRealType eigsTm;
		PsimagLite::diag(Tc, eigsTm, 'V');

		VectorComplexType expCoeff(mUsed, ComplexType(0));
		for (int k = 0; k < mUsed; ++k) {
			const ComplexType phase  = std::exp(ComplexType(0, -eigsTm[k] * dt));
			const ComplexType factor = std::conj(Tc(0, k)) * phase;
			for (int j = 0; j < mUsed; ++j)
				expCoeff[j] += factor * Tc(j, k);
		}

		VectorComplexType result(dim, ComplexType(0));
		for (int j = 0; j < mUsed; ++j) {
			const ComplexType c = norm0 * expCoeff[j];
			for (SizeType i = 0; i < dim; ++i)
				result[i] += c * Q[j][i];
		}
		return result;
	}

	// Exact exp(-i*H*dt)|psi> via dense diagonalization of csr (converted to a
	// full matrix). See krylovExpmvCSR's doc comment for why this is used
	// instead of the Lanczos recursion below a dim threshold.
	VectorComplexType
	denseExpmv(const CrsMatrixComplexType& csr, const VectorComplexType& psi, RealType dt) const
	{
		MatrixComplexType H = csr.toDense();
		VectorRealType    eigs;
		PsimagLite::diag(H, eigs, 'V');

		const SizeType    dim = psi.size();
		VectorComplexType coeff(dim, ComplexType(0));
		for (SizeType k = 0; k < dim; ++k) {
			ComplexType overlap(0);
			for (SizeType i = 0; i < dim; ++i)
				overlap += std::conj(H(i, k)) * psi[i];
			coeff[k] = overlap * std::exp(ComplexType(0, -eigs[k] * dt));
		}

		VectorComplexType result(dim, ComplexType(0));
		for (SizeType k = 0; k < dim; ++k)
			for (SizeType i = 0; i < dim; ++i)
				result[i] += H(i, k) * coeff[k];
		return result;
	}

	// One-shot equivalence check: assert max|applyHext(v) - CSR*v| < tol for both sectors.
	// Called once in solveLplus after buildHextCSR to validate the CSR construction.
	// Uses a non-zero artificial vMid to exercise the variable second-bath entries.
	void checkApplyHextEquivalence(ExtendedSector& sector) const
	{
		if (sector.upWordsNm1.empty() || sector.upWordsNp1.empty())
			return;

		// Artificial non-zero vMid so variable entries are non-trivially tested.
		std::vector<ComplexType> vMid(bathRank_);
		for (SizeType p = 0; p < bathRank_; ++p)
			vMid[p] = ComplexType(RealType(0.3) + RealType(0.1) * p,
			                      RealType(0.2) - RealType(0.05) * p);

		auto checkSector = [&](const std::string&           name,
		                       const std::vector<WordType>& upWords,
		                       const std::vector<WordType>& dnWords,
		                       SizeType                     dim1,
		                       CrsMatrixComplexType&        csr,
		                       const std::vector<VarEntry>& varEntries)
		{
			const SizeType dim = upWords.size() * dnWords.size();
			// Normalised unit vector as test input.
			VectorComplexType testV(
			    dim, ComplexType(RealType(1) / std::sqrt(RealType(dim))));

			updateCSR(csr, varEntries, vMid);

			VectorComplexType hvRef(dim, ComplexType(0));
			applyHext(testV, hvRef, vMid, upWords, dnWords, dim1);

			VectorComplexType hvCsr(dim, ComplexType(0));
			csr.matrixVectorProduct(hvCsr, testV);

			RealType maxDiff = 0;
			for (SizeType i = 0; i < dim; ++i)
				maxDiff = std::max(maxDiff, std::abs(hvRef[i] - hvCsr[i]));

			std::cout << "  checkApplyHextEquivalence [" << name << "]: dim=" << dim
			          << "  max|applyHext - CSR*v| = " << maxDiff << "\n";
			assert(maxDiff < RealType(1e-12)
			       && "applyHextCSR disagrees with applyHext — CSR build error");
		};

		checkSector("N-1",
		            sector.upWordsNm1,
		            sector.dnWordsNm1,
		            sector.dim1Nm1,
		            sector.csrNm1,
		            sector.varNm1);
		checkSector("N+1",
		            sector.upWordsNp1,
		            sector.dnWordsNp1,
		            sector.dim1Np1,
		            sector.csrNp1,
		            sector.varNp1);
		checkSector(
		    "N", sector.upWordsN, sector.dnWordsN, sector.dim1N, sector.csrN, sector.varN);
	}

	// ========== Fock-space Hamiltonian matrix-vector product ==========

	// Apply the post-quench extended Hamiltonian H_ext(V^+) to state v → Hv.
	// H = U_f n↑n↓|imp + Σ_i potPost_[i]*(n↑+n↓)|_i
	//   + Σ_{α,σ} V_α (c†_0 c_{α+1} + h.c.)             [first bath, real]
	//   + Σ_{p,σ} V^+_p (c†_0 c_{ep} + h.c.)            [second bath empty sites]
	//   + Σ_{p,σ} V^+_p (c†_0 c_{op} + h.c.)            [second bath occupied sites]
	void applyHext(const VectorComplexType&        v,
	               VectorComplexType&              hv,
	               const std::vector<ComplexType>& gbekHop,
	               const std::vector<WordType>&    upWords,
	               const std::vector<WordType>&    dnWords,
	               SizeType                        dim1) const
	{
		const SizeType dim    = v.size();
		const SizeType spinUp = LanczosPlusPlus::LanczosGlobals::SPIN_UP;
		const SizeType dim2   = (dim1 > 0) ? (dim / dim1) : 0;

		std::fill(hv.begin(), hv.end(), ComplexType(0));

		for (SizeType m = 0; m < dim; ++m) {
			if (std::abs(v[m]) < RealType(1e-15))
				continue;

			const SizeType    y    = m / dim1;
			const SizeType    x    = m % dim1;
			const WordType    up_m = upWords[x];
			const WordType    dn_m = dnWords[y];
			const ComplexType vm   = v[m];

			// --- Diagonal terms ---
			RealType diag = 0;
			if ((up_m & WordType(1)) && (dn_m & WordType(1)))
				diag += params_.uFinal; // Hubbard U at impurity (site 0)
			for (SizeType i = 0; i < nsites_ext_; ++i) {
				if ((up_m >> i) & WordType(1))
					diag += potPost_[i];
				if ((dn_m >> i) & WordType(1))
					diag += potPost_[i];
			}
			hv[m] += ComplexType(diag) * vm;

			// --- Hopping terms for each spin ---
			for (SizeType spin = 0; spin < 2; ++spin) {
				const WordType w       = (spin == spinUp) ? up_m : dn_m;
				const WordType w_other = (spin == spinUp) ? dn_m : up_m;

				// Apply c†_to c_from (hop from → to) on spin word w
				auto applyHop = [&](SizeType from, SizeType to, ComplexType t)
				{
					if (!((w >> from) & WordType(1)))
						return; // source unoccupied
					if ((w >> to) & WordType(1))
						return; // target occupied

					const WordType new_w
					    = (w & ~(WordType(1) << from)) | (WordType(1) << to);
					const WordType new_up = (spin == spinUp) ? new_w : w_other;
					const WordType new_dn = (spin == spinUp) ? w_other : new_w;

					const SizeType m_new = lookupIndex(
					    new_up, new_dn, upWords, dnWords, dim1, dim2);
					if (m_new == std::numeric_limits<SizeType>::max())
						return;

					const int sign = computeJWSign(w, from, to);
					hv[m_new] += t * ComplexType(sign) * vm;
				};

				// First bath: site 0 ↔ site α+1 (real hoppings V_α)
				for (SizeType a = 0; a < nBath_; ++a) {
					const ComplexType Va(firstBathHop_[a]);
					applyHop(0, a + 1, Va);
					applyHop(a + 1, 0, Va); // real: conj = same
				}

				// Second bath: site 0 ↔ empty sites and occupied sites (complex
				// V^+_p)
				for (SizeType p = 0; p < bathRank_; ++p) {
					const SizeType    emptyS = 1 + nBath_ + p;
					const SizeType    occS   = 1 + nBath_ + bathRank_ + p;
					const ComplexType Vp     = gbekHop[p];
					applyHop(0, emptyS, Vp);
					applyHop(emptyS, 0, std::conj(Vp));
					// Occupied orbital has the OPPOSITE conjugation from the
					// empty orbital's -- see buildHextCSR's matching comment
					// (GBEK Eq. 49/52a).
					applyHop(0, occS, std::conj(Vp));
					applyHop(occS, 0, Vp);
				}
			}
		}
	}

	// ========== Fock-space index lookup ==========

	// Find the basis index for (up_word, dn_word) using pre-sorted arrays.
	// Returns std::numeric_limits<SizeType>::max() if not found (shouldn't happen for valid
	// hops).
	static SizeType lookupIndex(WordType                     up,
	                            WordType                     dn,
	                            const std::vector<WordType>& upWords,
	                            const std::vector<WordType>& dnWords,
	                            SizeType                     dim1,
	                            SizeType /*dim2*/)
	{
		const auto it = std::lower_bound(upWords.begin(), upWords.end(), up);
		if (it == upWords.end() || *it != up)
			return std::numeric_limits<SizeType>::max();
		const SizeType x = static_cast<SizeType>(it - upWords.begin());

		const auto jt = std::lower_bound(dnWords.begin(), dnWords.end(), dn);
		if (jt == dnWords.end() || *jt != dn)
			return std::numeric_limits<SizeType>::max();
		const SizeType y = static_cast<SizeType>(jt - dnWords.begin());

		return x + y * dim1;
	}

	// Jordan-Wigner sign for c†_{to} c_{from} applied to spin word w.
	// Assumes bit 'from' is set and bit 'to' is unset in w.
	static int computeJWSign(WordType w, SizeType from, SizeType to)
	{
		// P_from: # bits at positions 0..from-1 in w
		int parity = 0;
		for (SizeType k = 0; k < from; ++k)
			parity += static_cast<int>((w >> k) & WordType(1));
		// After destroying at 'from':
		const WordType w1 = w ^ (WordType(1) << from);
		// P_to: # bits at positions 0..to-1 in w1
		for (SizeType k = 0; k < to; ++k)
			parity += static_cast<int>((w1 >> k) & WordType(1));
		return (parity & 1) ? -1 : 1;
	}

	// ========== Green's function formulas ==========
	//
	// Both rows below are computed via ONE backward Krylov sweep starting
	// from the t_n-seeded state (bStates[n]/dStates[n]), rather than as
	// independent (n,j) pairs -- see propagateOneStep's doc comment for why
	// per-time reseeding (not a single t=0 seed) is required, and
	// cross_check_seed_scheme.py for the derivation/confirmation of the
	// formulas used here.

	// Row of G^<(t_n, t_j) for all j=0..n, one spin configuration.
	// Derivation: G^<(t,t') = i<c†(t')c(t)>, and in the Heisenberg picture
	// c(t)|ψ0> = U_{N-1}(0,t) b(t) with b(t) = c·φ_N(t) = bStates at that
	// time, giving G^<(t_n,t_j) = i<b(t_j)|U_{N-1}(t_j,t_n)|b(t_n)> for j<=n
	// -- i.e. propagate bStates[n] BACKWARD to each earlier t_j.
	std::vector<ComplexType> gLesserRowGBEKSector(ExtendedSector& sector, int n) const
	{
		std::vector<ComplexType> row(static_cast<SizeType>(n + 1));
		VectorComplexType        psi  = sector.bStates[static_cast<SizeType>(n)];
		row[static_cast<SizeType>(n)] = ComplexType(0, 1) * innerProduct(psi, psi);

		for (int k = n - 1; k >= 0; --k) {
			const std::vector<ComplexType> vMid = computeVMid(k + 1);
			updateCSR(sector.csrNm1, sector.varNm1, vMid);
			psi = krylovExpmvCSR(psi, sector.csrNm1, -params_.dt);
			row[static_cast<SizeType>(k)] = ComplexType(0, 1)
			    * innerProduct(sector.bStates[static_cast<SizeType>(k)], psi);
		}
		return row;
	}

	// Row of G^>(t_n, t_j) for all j=0..n, one spin configuration.
	// Derivation: G^>(t,t') = -i<c(t)c†(t')>, giving
	// G^>(t_n,t_j) = -i<d(t_n)|U_{N+1}(t_n,t_j)|d(t_j)>
	//              = -i * conj( <d(t_j)|U_{N+1}(t_j,t_n)|d(t_n)> )   for j<=n,
	// so the SAME backward sweep (from dStates[n], through the N+1-sector
	// Hamiltonian) applies -- only the diagonal prefactor and the final
	// conjugation differ from gLesserRowGBEKSector.
	std::vector<ComplexType> gGreaterRowGBEKSector(ExtendedSector& sector, int n) const
	{
		std::vector<ComplexType> row(static_cast<SizeType>(n + 1));
		VectorComplexType        psi  = sector.dStates[static_cast<SizeType>(n)];
		row[static_cast<SizeType>(n)] = ComplexType(0, -1) * innerProduct(psi, psi);

		for (int k = n - 1; k >= 0; --k) {
			const std::vector<ComplexType> vMid = computeVMid(k + 1);
			updateCSR(sector.csrNp1, sector.varNp1, vMid);
			psi = krylovExpmvCSR(psi, sector.csrNp1, -params_.dt);
			const ComplexType inner
			    = innerProduct(sector.dStates[static_cast<SizeType>(k)], psi);
			row[static_cast<SizeType>(k)] = ComplexType(0, -1) * std::conj(inner);
		}
		return row;
	}

	// ========== Double occupation and energy observables (paper Figs. 9-10) ==========
	//
	// <d(t)> = <n_{imp,up}(t) n_{imp,dn}(t)>, a diagonal (time-local) operator
	// expectation value on PhiNHist[n] -- no Hamiltonian needed at all.
	SizeType
	basisIndexToXY(const ExtendedSector& sector, SizeType m, SizeType& x, SizeType& y) const
	{
		x = m % sector.dim1N;
		y = m / sector.dim1N;
		return m;
	}

	RealType computeDoubleOccupationGBEKSector(const ExtendedSector& sector, int n) const
	{
		const VectorComplexType& psi  = sector.PhiNHist[static_cast<SizeType>(n)];
		RealType                 docc = 0;
		for (SizeType m = 0; m < psi.size(); ++m) {
			SizeType x = 0, y = 0;
			basisIndexToXY(sector, m, x, y);
			const WordType up_m = sector.upWordsN[x];
			const WordType dn_m = sector.dnWordsN[y];
			if ((up_m & WordType(1)) && (dn_m & WordType(1)))
				docc += std::norm(psi[m]);
		}
		return docc;
	}

	// <Ekin(t_n)> = <PhiNHist[n]|H_hop(t_mid)|PhiNHist[n]>, where H_hop is
	// csrN with its fixed onsite (U+potential) diagonal removed. PRECONDITION:
	// sector.csrN must already hold H(t_mid) for the interval ending at step
	// n -- true whenever this is called right after ensurePropagated(n) (as
	// computeGimp does), since propagateOneStep(n) is the last thing to have
	// called updateCSR(sector.csrN, ...) for step n (for n=0, csrN instead
	// holds its buildHextCSR construction-time state, i.e. V=0 -- correct,
	// since Ekin(0)=0 exactly before any hopping switches on). Deliberately
	// does NOT recompute/reload vMid itself: doing so would mean calling
	// computeVMid(n), which reads decomp_->Vplus(...) -- wrong whenever csrN
	// was instead driven directly (bypassing NeqBathDecomposition entirely),
	// as the hand-driven Catch2 tests do.
	//
	// The 0.5 factor: the bare <H_hop> expectation value comes out at
	// exactly 2x the physical kinetic energy. Each bath pair couples to
	// the impurity through TWO auxiliary sites (an "occ" partner carrying
	// Lambda^< and an "empty" partner carrying Lambda^>, per Fig. 5/Eq.
	// 56-63), each with the same coupling magnitude |V_p(t)| -- correct
	// for reproducing Lambda_+ itself (confirmed independently: d(t)
	// matches the paper's Fig. 10 exactly) but it double-counts the single
	// physical hybridization channel when read directly off <H_hop>.
	// Caught in the Python reference by cross-checking against the exact
	// conservation law <Etot(t)> = <Etot(0)> = -U/4 for t > t_q, which the
	// raw (uncorrected) value badly violated -- see
	// cincuenta/TestSuite/gbek_reference/gbek_selfconsistency.py
	// (compute_energy_observables) and project memory
	// project_gbek_ekin_factor2_bugfix. The conservation-law regression
	// test in test_ImpuritySolverNeqGBEK.cpp guards this factor directly.
	RealType computeKineticEnergyGBEKSector(const ExtendedSector& sector, int n) const
	{
		assert(n >= 0 && n < static_cast<int>(sector.PhiNHist.size()));
		const SizeType nn = static_cast<SizeType>(n);

		const VectorComplexType& psi = sector.PhiNHist[nn];
		VectorComplexType        Hpsi;
		sparseMatVec(sector.csrN, psi, Hpsi);
		const ComplexType raw = innerProduct(psi, Hpsi);

		RealType onsite = 0;
		for (SizeType m = 0; m < psi.size(); ++m)
			onsite += sector.diagFixed[m] * std::norm(psi[m]);

		const ComplexType diff = raw - ComplexType(onsite);
		assert(std::abs(std::imag(diff)) < RealType(1e-6)
		       && "computeKineticEnergyGBEKSector: <H_hop> should be real "
		          "(expect ~0 imaginary part for a physical energy)");
		return RealType(0.5) * std::real(diff);
	}

	// Gα/Gβ average (GBEK Eq. 70): G_up = (1/2)(Gα_up + Gβ_up). When
	// nup_==ndown_, system beta is identical to system alpha and was never
	// built, so this reduces to just system alpha's row.
	std::vector<ComplexType> gLesserRowGBEK(int n) const
	{
		std::vector<ComplexType> rowA = gLesserRowGBEKSector(sectorAlpha_, n);
		if (sameConfig_)
			return rowA;
		const std::vector<ComplexType> rowB = gLesserRowGBEKSector(sectorBeta_, n);
		for (SizeType k = 0; k < rowA.size(); ++k)
			rowA[k] = RealType(0.5) * (rowA[k] + rowB[k]);
		return rowA;
	}

	std::vector<ComplexType> gGreaterRowGBEK(int n) const
	{
		std::vector<ComplexType> rowA = gGreaterRowGBEKSector(sectorAlpha_, n);
		if (sameConfig_)
			return rowA;
		const std::vector<ComplexType> rowB = gGreaterRowGBEKSector(sectorBeta_, n);
		for (SizeType k = 0; k < rowA.size(); ++k)
			rowA[k] = RealType(0.5) * (rowA[k] + rowB[k]);
		return rowA;
	}

	// dα/dβ average, same structure as gLesserRowGBEK/gGreaterRowGBEK above.
	RealType computeDoccGBEK(int n) const
	{
		const RealType dA = computeDoubleOccupationGBEKSector(sectorAlpha_, n);
		if (sameConfig_)
			return dA;
		const RealType dB = computeDoubleOccupationGBEKSector(sectorBeta_, n);
		return RealType(0.5) * (dA + dB);
	}

	RealType computeKineticEnergyGBEK(int n) const
	{
		const RealType eA = computeKineticEnergyGBEKSector(sectorAlpha_, n);
		if (sameConfig_)
			return eA;
		const RealType eB = computeKineticEnergyGBEKSector(sectorBeta_, n);
		return RealType(0.5) * (eA + eB);
	}

	// ========== Vector helpers ==========

	static void axpy(VectorComplexType& w, ComplexType c, const VectorComplexType& v)
	{
		for (SizeType i = 0; i < w.size(); ++i)
			w[i] += c * v[i];
	}

	static void scaleVec(VectorComplexType& v, ComplexType c)
	{
		for (SizeType i = 0; i < v.size(); ++i)
			v[i] *= c;
	}

	static ComplexType innerProduct(const VectorComplexType& a, const VectorComplexType& b)
	{
		ComplexType s(0);
		for (SizeType i = 0; i < a.size(); ++i)
			s += std::conj(a[i]) * b[i];
		return s;
	}

	// ========== LanczosPlusPlus infrastructure (mirrors ExactDiag) ==========

	// Minimal LanczosPlusPlus input string for a star-geometry Anderson model.
	static std::string buildLanczosInput(RealType              U,
	                                     SizeType              nup,
	                                     SizeType              ndown,
	                                     const VectorRealType& hoppings,
	                                     const VectorRealType& potV,
	                                     SizeType              nsites)
	{
		std::string uStr = "[" + ttos(U);
		for (SizeType i = 1; i < nsites; ++i)
			uStr += ", 0.";
		uStr += "]";

		std::string connStr = "[";
		for (SizeType i = 0; i < hoppings.size(); ++i) {
			if (i > 0)
				connStr += ",";
			connStr += ttos(hoppings[i]);
		}
		connStr += "]";

		std::string potStr = "[";
		for (SizeType i = 0; i < nsites; ++i) {
			if (i > 0)
				potStr += ",";
			potStr += ttos(potV[i]);
		}
		potStr += ",";
		for (SizeType i = 0; i < nsites; ++i) {
			if (i > 0)
				potStr += ",";
			potStr += ttos(potV[i]);
		}
		potStr += "]";

		std::string s = "##Ainur1.0\n\n";
		s += "TotalNumberOfSites=" + ttos(nsites) + ";\n";
		s += "NumberOfTerms=1;\n";
		s += "DegreesOfFreedom=1;\n";
		s += "GeometryKind=star;\n";
		s += "GeometryOptions=none;\n";
		s += "hubbardU=" + uStr + ";\n";
		s += "Model=HubbardOneBand;\n";
		s += "SolverOptions=twositedmrg,geometryallinsystem,hd5dontprint;\n";
		s += "Version=templateForDMFT;\n";
		s += "OutputFile=neqGbekDummy;\n";
		s += "InfiniteLoopKeptStates=1;\n";
		s += "FiniteLoops=0 0 0;\n";
		s += "TargetElectronsUp=" + ttos(nup) + ";\n";
		s += "TargetElectronsDown=" + ttos(ndown) + ";\n";
		s += "dir0:Connectors=" + connStr + ";\n";
		s += "potentialV=" + potStr + ";\n";
		return s;
	}

	// Full diagonalisation of a sector of the Hamiltonian.
	static void diagWithBasis(const ModelBaseType& model,
	                          const BasisBaseType& basis,
	                          const GeometryType&  geom,
	                          VectorRealType&      eigs,
	                          MatrixType&          eigvecs)
	{
		DefaultSymmetryType       rs(basis, geom, "");
		InternalProductStoredType ham(model, basis, rs);
		const SizeType            dim = basis.size();
		eigs.resize(dim);
		eigvecs.resize(dim, dim);
		ham.fullDiag(eigs, eigvecs);
	}

	// Basis index m such that basis(m,SPIN_UP)==occupiedSitesUp and
	// basis(m,SPIN_DOWN)==occupiedSitesDown, or SizeType max if not found.
	// Linear scan over the (already sector-restricted) basis -- dim is at
	// most a few hundred thousand here and this runs once per buildSector
	// call, so an O(dim) scan is negligible; avoids depending on the
	// sorted-array lookup machinery (buildFockLookup/lookupIndex) that
	// isn't built yet at this point in buildSector.
	static SizeType findBasisIndex(const BasisBaseType& basis,
	                               WordType             occupiedSitesUp,
	                               WordType             occupiedSitesDown)
	{
		const SizeType spinUp = LanczosPlusPlus::LanczosGlobals::SPIN_UP;
		const SizeType spinDn = LanczosPlusPlus::LanczosGlobals::SPIN_DOWN;
		for (SizeType m = 0; m < basis.size(); ++m)
			if (basis(m, spinUp) == occupiedSitesUp
			    && basis(m, spinDn) == occupiedSitesDown)
				return m;
		return std::numeric_limits<SizeType>::max();
	}

	// Lanczos ground state for sectors too large for full diagonalisation
	// (nsites_ext_>8). Returns a dim×1 matrix whose column 0 is the
	// normalised GS eigenvector.
	//
	// Seeded with the KNOWN atomic-limit product state the pre-quench
	// bigEps potentials (see buildSector's potPre construction) are
	// designed to make overwhelmingly the ground state: the impurity
	// carries impurityUp/impurityDown electrons (0 or 1 each), each of the
	// L "occ" second-bath sites (index 1+nBath+L+p) is doubly occupied,
	// every other site (including all nBath first-bath sites, whose true
	// equilibrium occupation isn't known in closed form here) is left
	// empty in the seed -- a neutral default for any first-bath sites,
	// though the only case actually exercised so far (NeqAtomicLimit=1) has
	// nBath=0, for which the pre-quench Hamiltonian has NO hopping at all,
	// so this seed is not just a good guess but the EXACT eigenstate.
	//
	// This replaces an earlier version that seeded Lanczos with a naive
	// uniform vector (1/sqrt(dim) on every basis state), confirmed to
	// converge to the WRONG extremal state for L>=4 (reported ground
	// energy -2807.96 vs. the true -4000 for L=4's fully-polarized product
	// state) -- see project memory project_gbek_cpp_lanczosGS_bug for the
	// full diagnosis. A seed with substantial overlap with the true ground
	// state, rather than an unbiased guess, is the fix.
	static MatrixType lanczosGS(const ModelBaseType& model,
	                            const GeometryType&  geom,
	                            SizeType             nBath,
	                            SizeType             L,
	                            SizeType             impurityUp,
	                            SizeType             impurityDown)
	{
		const BasisBaseType&      basis = model.basis();
		const SizeType            dim   = basis.size();
		DefaultSymmetryType       rs(basis, geom, "");
		InternalProductStoredType ham(model, basis, rs);

		PsimagLite::ParametersForSolver<RealType> lparams;
		lparams.steps      = static_cast<SizeType>(std::min(dim, SizeType(200)));
		lparams.tolerance  = 1e-12;
		lparams.lotaMemory = true;

		LanczosSolverForGSType lanczos(ham, lparams);

		// occupiedSitesUp/occupiedSitesDown: the set of sites occupied by an
		// up/down electron respectively (packed one bit per site), same
		// representation and same site layout as buildHextCSR (site 0 =
		// impurity; 1..nBath = first bath; nBath+1..nBath+L = "empty"
		// second-bath sites; nBath+L+1..nBath+2L = "occ" second-bath sites --
		// occS below matches that function's own occS computation exactly).
		// Marking occS occupied in BOTH sets means "this occ site holds one
		// up AND one down electron", i.e. doubly occupied.
		WordType occupiedSitesUp = 0, occupiedSitesDown = 0;
		if (impurityUp == 1)
			occupiedSitesUp |= (WordType(1) << 0);
		if (impurityDown == 1)
			occupiedSitesDown |= (WordType(1) << 0);
		for (SizeType p = 0; p < L; ++p) {
			const SizeType occS = 1 + nBath + L + p;
			occupiedSitesUp |= (WordType(1) << occS);
			occupiedSitesDown |= (WordType(1) << occS);
		}
		const SizeType seedIndex
		    = findBasisIndex(basis, occupiedSitesUp, occupiedSitesDown);

		typename LanczosSolverForGSType::VectorType initVec(dim, ComplexOrRealType(0));
		if (seedIndex != std::numeric_limits<SizeType>::max()) {
			initVec[seedIndex] = ComplexOrRealType(1);
		} else {
			assert(false
			       && "lanczosGS: target atomic-limit product state not found in "
			          "basis -- falling back to a uniform seed, which is known to "
			          "converge to the wrong state for large nsites_ext_");
			for (SizeType i = 0; i < dim; ++i)
				initVec[i] = ComplexOrRealType(1.0 / std::sqrt(RealType(dim)));
		}

		RealType                                    gsEnergy = 0;
		typename LanczosSolverForGSType::VectorType gsVec(dim, ComplexOrRealType(0));

		lanczos.computeOneState(gsEnergy, gsVec, initVec, 0);

		MatrixType result(dim, 1);
		for (SizeType i = 0; i < dim; ++i)
			result(i, 0) = gsVec[i];
		return result;
	}

	// ========== Member variables ==========

	SizeType                    bathRank_;
	const ParamsNeqType&        params_;
	ExactDiagType               exactDiag_;
	std::unique_ptr<DecompType> decomp_;

	// Electron counts (read from io for L>0)
	SizeType nup_, ndown_;

	// Extended system geometry
	SizeType nBath_; // # first-bath sites
	SizeType nsites_ext_; // total sites = nBath+1+2L

	// First bath parameters (extracted from bathParams in solveLplus)
	VectorRealType firstBathHop_;
	VectorRealType firstBathEps_;

	// Post-quench on-site potentials for the extended system (shared by both
	// spin configurations -- U_final and first-bath epsilons don't depend on
	// which spin the impurity's extra electron carries).
	VectorRealType potPost_;

	// Gα (nup_ext=nup+L, ndown_ext=ndown+L) and Gβ (swapped) extended-Fock
	// state (GBEK Eq. 70). sameConfig_ is true when nup_==ndown_, in which
	// case sectorBeta_ is never built and is identical to sectorAlpha_.
	mutable ExtendedSector sectorAlpha_, sectorBeta_;
	bool                   sameConfig_ = false;

	// Double occupation and kinetic energy history, filled progressively by
	// computeGimp(n) alongside the Green's functions (one entry per time
	// step actually visited); see dumpDoccAndEnergy for the derived
	// Eint/Etot columns and paper Figs. 9-10.
	mutable std::vector<RealType> doccHistory_, ekinHistory_;

	friend struct GBEKTestAccessor;
};

} // namespace Dmft
#endif // IMPURITYSOLVER_NEQ_GBEK_H
