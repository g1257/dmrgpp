#pragma once
// Shared test utilities for cincuenta unit tests.
//
// GBEKTestAccessor — named friend accessor for private members of
//   ImpuritySolverNeqGBEK<std::complex<double>>.
//
// assertHermitian, buildDenseFromApplyHext, buildDenseFromCSR, vecNorm —
//   matrix/vector check utilities reusable across GBEK and ExactDiag tests.

// ImpuritySolverNeqGBEK.h transitively pulls in CrsMatrix.h, Matrix.h, etc.
#include "ImpuritySolverNeqGBEK.h"
#include <catch2/catch_approx.hpp>
#include <catch2/catch_test_macros.hpp>

namespace Dmft {

// ── GBEKTestAccessor ────────────────────────────────────────────────────────
// Named friend accessor for ImpuritySolverNeqGBEK<std::complex<double>>.
// The class declares: friend struct GBEKTestAccessor;

struct GBEKTestAccessor {
	using SolverType    = ImpuritySolverNeqGBEK<std::complex<double>>;
	using RealType      = double;
	using ComplexType   = std::complex<double>;
	using CrsMatrixType = PsimagLite::CrsMatrix<ComplexType>;
	using VarEntry      = typename SolverType::VarEntry;
	using VectorComplex = typename PsimagLite::Vector<ComplexType>::Type;
	using WordType      = LanczosPlusPlus::LanczosGlobals::WordType;

	// ── Const data-member accessors ──────────────────────────────────────
	// All route through sectorAlpha_ (system alpha; identical to the old
	// single-configuration fields when nup_==ndown_, which is what the
	// existing tests exercise).
	static const CrsMatrixType& csrNm1(const SolverType& s) { return s.sectorAlpha_.csrNm1; }
	static const CrsMatrixType& csrNp1(const SolverType& s) { return s.sectorAlpha_.csrNp1; }
	// mutable CSR exposed for updateCSR tests
	static CrsMatrixType& csrNm1Mut(const SolverType& s) { return s.sectorAlpha_.csrNm1; }
	static CrsMatrixType& csrNp1Mut(const SolverType& s) { return s.sectorAlpha_.csrNp1; }

	static const std::vector<VarEntry>& varNm1(const SolverType& s)
	{
		return s.sectorAlpha_.varNm1;
	}
	static const std::vector<VarEntry>& varNp1(const SolverType& s)
	{
		return s.sectorAlpha_.varNp1;
	}

	static const std::vector<WordType>& upWordsNm1(const SolverType& s)
	{
		return s.sectorAlpha_.upWordsNm1;
	}
	static const std::vector<WordType>& dnWordsNm1(const SolverType& s)
	{
		return s.sectorAlpha_.dnWordsNm1;
	}
	static SizeType dim1Nm1(const SolverType& s) { return s.sectorAlpha_.dim1Nm1; }

	static const std::vector<WordType>& upWordsNp1(const SolverType& s)
	{
		return s.sectorAlpha_.upWordsNp1;
	}
	static const std::vector<WordType>& dnWordsNp1(const SolverType& s)
	{
		return s.sectorAlpha_.dnWordsNp1;
	}
	static SizeType dim1Np1(const SolverType& s) { return s.sectorAlpha_.dim1Np1; }

	static const std::vector<VectorComplex>& bStates(const SolverType& s)
	{
		return s.sectorAlpha_.bStates;
	}
	static const std::vector<VectorComplex>& phiNHist(const SolverType& s)
	{
		return s.sectorAlpha_.PhiNHist;
	}
	static SizeType bathRank(const SolverType& s) { return s.bathRank_; }

	// N-sector Hamiltonian CSR + the c_{imp,up} operator (N -> N-1 sector),
	// needed to reproduce the reseed-at-every-step construction directly
	// (see the seed-scheme regression test).
	static const CrsMatrixType& csrN(const SolverType& s) { return s.sectorAlpha_.csrN; }
	static CrsMatrixType&       csrNMut(const SolverType& s) { return s.sectorAlpha_.csrN; }
	static const std::vector<VarEntry>& varN(const SolverType& s)
	{
		return s.sectorAlpha_.varN;
	}
	static const CrsMatrixType& cUpNm1(const SolverType& s) { return s.sectorAlpha_.cUpNm1; }
	static const std::vector<WordType>& upWordsN(const SolverType& s)
	{
		return s.sectorAlpha_.upWordsN;
	}
	static const std::vector<WordType>& dnWordsN(const SolverType& s)
	{
		return s.sectorAlpha_.dnWordsN;
	}
	static SizeType dim1N(const SolverType& s) { return s.sectorAlpha_.dim1N; }

	// Mutable access to PhiNHist, for tests that hand-drive a chosen V(t)
	// sequence directly (bypassing NeqBathDecomposition/self-consistency)
	// and need to stash the resulting states before calling the docc/energy
	// observable methods below -- same style as csrNMut.
	static std::vector<VectorComplex>& phiNHistMut(const SolverType& s)
	{
		return s.sectorAlpha_.PhiNHist;
	}
	static const std::vector<RealType>& diagFixed(const SolverType& s)
	{
		return s.sectorAlpha_.diagFixed;
	}

	// docc(t_n)/Ekin(t_n), alpha/beta-averaged (paper Figs. 9-10).
	static RealType computeDoccGBEK(const SolverType& s, int n) { return s.computeDoccGBEK(n); }
	static RealType computeKineticEnergyGBEK(const SolverType& s, int n)
	{
		return s.computeKineticEnergyGBEK(n);
	}

	// Single-sector (system alpha only) versions, for tests that hand-drive
	// only sectorAlpha_'s PhiNHist directly (see phiNHistMut) and don't want
	// to also drive sectorBeta_ just to exercise the docc/Ekin arithmetic
	// itself. Physically valid on their own (conservation of <psi|H(t)|psi>
	// under unitary evolution holds for any single trajectory), just not the
	// paper's Eq.-70-symmetrized quantity.
	static RealType computeDoubleOccupationGBEKSector(const SolverType& s, int n)
	{
		return s.computeDoubleOccupationGBEKSector(s.sectorAlpha_, n);
	}
	static RealType computeKineticEnergyGBEKSector(const SolverType& s, int n)
	{
		return s.computeKineticEnergyGBEKSector(s.sectorAlpha_, n);
	}
	static const std::vector<RealType>& doccHistory(const SolverType& s)
	{
		return s.doccHistory_;
	}
	static const std::vector<RealType>& ekinHistory(const SolverType& s)
	{
		return s.ekinHistory_;
	}

	// ── Method delegators ────────────────────────────────────────────────
	static void applyHext(const SolverType&               s,
	                      const VectorComplex&            v,
	                      VectorComplex&                  hv,
	                      const std::vector<ComplexType>& gbekHop,
	                      const std::vector<WordType>&    upWords,
	                      const std::vector<WordType>&    dnWords,
	                      SizeType                        dim1)
	{
		s.applyHext(v, hv, gbekHop, upWords, dnWords, dim1);
	}

	static VectorComplex krylovExpmv(const SolverType&               s,
	                                 const VectorComplex&            psi,
	                                 const std::vector<ComplexType>& vMid,
	                                 const std::vector<WordType>&    upWords,
	                                 const std::vector<WordType>&    dnWords,
	                                 SizeType                        dim1,
	                                 RealType                        dt)
	{
		return s.krylovExpmv(psi, vMid, upWords, dnWords, dim1, dt);
	}

	static VectorComplex krylovExpmvCSR(const SolverType&    s,
	                                    const VectorComplex& psi,
	                                    const CrsMatrixType& csr,
	                                    RealType             dt)
	{
		return s.krylovExpmvCSR(psi, csr, dt);
	}

	static void updateCSR(const SolverType&               s,
	                      CrsMatrixType&                  csr,
	                      const std::vector<VarEntry>&    varEntries,
	                      const std::vector<ComplexType>& vMid)
	{
		s.updateCSR(csr, varEntries, vMid);
	}

	// Rectangular-safe sparse mat-vec (x = A*y), for cUpNm1 (N -> N-1 sector,
	// not square).
	static VectorComplex
	sparseMatVec(const SolverType& s, const CrsMatrixType& A, const VectorComplex& y)
	{
		VectorComplex x;
		s.sparseMatVec(A, y, x);
		return x;
	}
};

} // namespace Dmft

// ── Free utility functions ───────────────────────────────────────────────────

using GBEKUtils_Complex   = std::complex<double>;
using GBEKUtils_Real      = double;
using GBEKUtils_Matrix    = PsimagLite::Matrix<GBEKUtils_Complex>;
using GBEKUtils_CrsMatrix = PsimagLite::CrsMatrix<GBEKUtils_Complex>;
using GBEKUtils_Vector    = typename PsimagLite::Vector<GBEKUtils_Complex>::Type;

// Assert max|M - M†| < tol; REQUIRE square, CHECK max error.
inline void assertHermitian(const GBEKUtils_Matrix& M, GBEKUtils_Real tol = 1e-12)
{
	REQUIRE(M.n_row() == M.n_col());
	const SizeType n      = M.n_row();
	GBEKUtils_Real maxErr = 0;
	SizeType       wi = 0, wj = 0;
	for (SizeType i = 0; i < n; ++i) {
		for (SizeType j = 0; j < n; ++j) {
			const GBEKUtils_Real e = std::abs(M(i, j) - std::conj(M(j, i)));
			if (e > maxErr) {
				maxErr = e;
				wi     = i;
				wj     = j;
			}
		}
	}
	INFO("max|H[i,j]-conj(H[j,i])| = " << maxErr << " at (" << wi << "," << wj << ")");
	CHECK(maxErr < tol);
}

// Build a dense Hilbert-space matrix from applyHext: H[:,j] = applyHext(e_j).
template <typename SolverType>
GBEKUtils_Matrix
buildDenseFromApplyHext(const SolverType&                                             solver,
                        const std::vector<GBEKUtils_Complex>&                         gbekHop,
                        const std::vector<LanczosPlusPlus::LanczosGlobals::WordType>& upWords,
                        const std::vector<LanczosPlusPlus::LanczosGlobals::WordType>& dnWords,
                        SizeType                                                      dim1)
{
	const SizeType   dim = upWords.size() * dnWords.size();
	GBEKUtils_Matrix H(dim, dim, GBEKUtils_Complex(0));
	GBEKUtils_Vector ej(dim, GBEKUtils_Complex(0));
	GBEKUtils_Vector hv(dim, GBEKUtils_Complex(0));

	for (SizeType j = 0; j < dim; ++j) {
		ej[j] = GBEKUtils_Complex(1);
		std::fill(hv.begin(), hv.end(), GBEKUtils_Complex(0));
		Dmft::GBEKTestAccessor::applyHext(solver, ej, hv, gbekHop, upWords, dnWords, dim1);
		for (SizeType i = 0; i < dim; ++i)
			H(i, j) = hv[i];
		ej[j] = GBEKUtils_Complex(0);
	}
	return H;
}

// Convert a CrsMatrix to a dense PsimagLite Matrix.
inline GBEKUtils_Matrix buildDenseFromCSR(const GBEKUtils_CrsMatrix& csr)
{
	GBEKUtils_Matrix M(csr.rows(), csr.cols(), GBEKUtils_Complex(0));
	crsMatrixToFullMatrix(M, csr);
	return M;
}

// Euclidean norm of a complex vector.
inline GBEKUtils_Real vecNorm(const GBEKUtils_Vector& v)
{
	GBEKUtils_Real s = 0;
	for (const auto& x : v)
		s += std::norm(x);
	return std::sqrt(s);
}
