// #include "BatchedGemmKokkos.h"
#include "BatchedGemmCpu.h"
#include "BatchedGemmPluginSc.h"
#include <Kokkos_Core.hpp>
#include <PsimagLite/CrsMatrix.h>
#include <PsimagLite/Matrix.h>
#include <PsimagLite/Vector.h>
#include <catch2/catch_approx.hpp>
#include <catch2/catch_test_macros.hpp>
#include <random>

using SizeType = unsigned long;

// Minimal fake types to instantiate BatchedGemmKokkos for a 1-patch, 1-operator case
struct FakeParams {
	struct Opt {
		bool isSet(const std::string& s) const { return s == "BatchedGemm"; }
	} options;
};

struct FakeLR {
	std::vector<int> parts; // partition points, length = groups+1
	int              partition(size_t i) const { return parts[i]; }
	int              size() const { return parts.back(); }
};

struct FakeLrs {
	FakeLR        l, r;
	const FakeLR& left() const { return l; }
	const FakeLR& right() const { return r; }
};

// Matrix wrapper used by BatchedGemmKokkos via ArrayOfMatStructType.
struct FakeMatrixDenseOrSparse {
	using value_type = double;
	using VectorType = PsimagLite::Vector<value_type>::Type;
	using MatrixType = PsimagLite::Matrix<value_type>;
	using CrsType    = PsimagLite::CrsMatrix<value_type>;
	MatrixType mat;
	CrsType    crs;
	bool       zero = false;
	FakeMatrixDenseOrSparse()
	    : mat()
	    , crs(0, 0)
	    , zero(true)
	{ }
	FakeMatrixDenseOrSparse(const MatrixType& m)
	    : mat(m)
	    , crs(0, 0)
	    , zero(false)
	{ }
	bool              isZero() const { return zero; }
	bool              isDense() const { return true; }
	const MatrixType& dense() const { return mat; }
	const CrsType&    sparse() const { return crs; }
	int               rows() const { return mat.rows(); }
	int               cols() const { return mat.cols(); }
};

// ArrayOfMatStructType: callable (ip,jp) -> pointer to FakeMatrixDenseOrSparse
struct FakeArrayOfMatStructType {
	using MatrixDenseOrSparseType = FakeMatrixDenseOrSparse;

	std::vector<std::vector<FakeMatrixDenseOrSparse>> storage; // [ip][jp]
	FakeArrayOfMatStructType(size_t npatches = 0)
	{
		storage.resize(npatches);
		for (size_t i = 0; i < npatches; ++i)
			storage[i].resize(npatches);
	}
	const FakeMatrixDenseOrSparse* operator()(size_t ip, size_t jp) const
	{
		// If matrix is default-constructed, treat as zero -> return nullptr
		if (storage[ip][jp].isZero())
			return nullptr;
		return &storage[ip][jp];
	}
	FakeMatrixDenseOrSparse* operator()(size_t ip, size_t jp)
	{
		if (storage[ip][jp].isZero())
			return nullptr;
		return &storage[ip][jp];
	}
	void set(size_t ip, size_t jp, const FakeMatrixDenseOrSparse& m) { storage[ip][jp] = m; }
};

struct GenIjPatch {
	enum LeftOrRightEnumType
	{
		LEFT  = 0,
		RIGHT = 1
	};
	using BasisType = int;
};

struct FakeInitKron {
	enum WhatBasisEnum
	{
		OLD = 0,
		NEW = 1
	};

	using ArrayOfMatStructType = FakeArrayOfMatStructType;
	using GenIjPatchType       = GenIjPatch;
	using SparseMatrixType     = PsimagLite::CrsMatrix<double>;

	FakeParams               p;
	FakeArrayOfMatStructType xc0, yc0;
	FakeLrs                  lrs_; // used for both NEW/OLD in this fake
	std::vector<size_t>      patchLeft, patchRight;
	size_t                   npatches_  = 0;
	size_t                   noperator_ = 0;

	FakeInitKron(size_t npatches, size_t noperator)
	    : p()
	    , xc0(npatches)
	    , yc0(npatches)
	    , npatches_(npatches)
	    , noperator_(noperator)
	{
		// default partitions: single group covering full size; user will set lrs_.l.parts
		// and lrs_.r.parts
		patchLeft.resize(npatches_, 0);
		patchRight.resize(npatches_, 0);
	}

	const FakeParams& params() const { return p; }

	size_t numberOfPatches(int /*which*/) const { return npatches_; }
	size_t connections() const { return noperator_; }

	// Return patch group indices
	const std::vector<size_t>& patch(WhatBasisEnum /*which*/,
	                                 GenIjPatch::LeftOrRightEnumType side) const
	{
		if (side == GenIjPatch::LEFT)
			return patchLeft;
		return patchRight;
	}

	const FakeLrs& lrs(WhatBasisEnum /*which*/) const { return lrs_; }

	// offsetForPatches: cumulative offsets in vin/vout per patch
	SizeType offsetForPatches(WhatBasisEnum /*which*/, size_t ipatch) const
	{
		SizeType off = 0;
		for (size_t p = 0; p < ipatch; ++p) {
			size_t igroup = patchLeft[p];
			size_t jgroup = patchRight[p];
			int    L1     = lrs_.left().partition(igroup);
			int    L2     = lrs_.left().partition(igroup + 1);
			int    R1     = lrs_.right().partition(jgroup);
			int    R2     = lrs_.right().partition(jgroup + 1);
			off += static_cast<SizeType>((L2 - L1) * (R2 - R1));
		}
		return off;
	}

	// Accessors for matrices: expect xc(k) and yc(k) to provide (ip,jp) access
	const FakeArrayOfMatStructType& xc(size_t /*k*/) const { return xc0; }
	const FakeArrayOfMatStructType& yc(size_t /*k*/) const { return yc0; }

	void
	checks(const FakeMatrixDenseOrSparse&, const FakeMatrixDenseOrSparse&, size_t, size_t) const
	{ }
};

template <template <typename> class Impl> void runBatchedGemmTest(bool useKokkos)
{
	if (useKokkos)
		Kokkos::initialize();

	// Build fake initKron for 1 patch, 1 operator
	FakeInitKron fk(1, 1);

	// Define partitions: left size = 2, right size = 3
	fk.lrs_.l.parts = { 0, 2 }; // partition points: 0..2 -> size = 2
	fk.lrs_.r.parts = { 0, 3 }; // size = 3

	fk.patchLeft[0]  = 0; // trivial group index
	fk.patchRight[0] = 0; // trivial group index for right

	// Random generator with fixed seed for reproducibility
	std::mt19937_64                        rng(123456789);
	std::uniform_real_distribution<double> dist(-1.0, 1.0);

	// Create random matrices A (2x2) and B (3x3)
	PsimagLite::Matrix<double> A(2, 2);
	for (int i = 0; i < 2; ++i)
		for (int j = 0; j < 2; ++j)
			A(i, j) = dist(rng);

	PsimagLite::Matrix<double> B(3, 3);
	for (int i = 0; i < 3; ++i)
		for (int j = 0; j < 3; ++j)
			B(i, j) = dist(rng);

	fk.xc0.set(0, 0, FakeMatrixDenseOrSparse(A));
	fk.yc0.set(0, 0, FakeMatrixDenseOrSparse(B));

	// Instantiate implementation inside a scope so it is destroyed before Kokkos::finalize()
	{
		Impl<FakeInitKron> bg(fk);

		// Prepare vin as column-major 3x2 with random entries
		using VectorType = PsimagLite::Vector<double>::Type;
		VectorType vin(6), vout(6);
		for (size_t i = 0; i < 6; ++i)
			vin[i] = dist(rng);
		for (size_t i = 0; i < 6; ++i)
			vout[i] = 0.0;

		// Compute expected result on host
		PsimagLite::Matrix<double> X(3, 2);
		for (size_t col = 0; col < 2; ++col)
			for (size_t row = 0; row < 3; ++row)
				X(row, col) = vin[row + col * 3];

		PsimagLite::Matrix<double> BX(3, 2);
		for (size_t i = 0; i < 3; ++i)
			for (size_t j = 0; j < 2; ++j) {
				double s = 0.0;
				for (size_t k = 0; k < 3; ++k)
					s += B(i, k) * X(k, j);
				BX(i, j) = s;
			}

		PsimagLite::Matrix<double> Y(3, 2);
		for (size_t i = 0; i < 3; ++i)
			for (size_t j = 0; j < 2; ++j) {
				double s = 0.0;
				for (size_t k = 0; k < 2; ++k)
					s += BX(i, k) * A(j, k);
				Y(i, j) = s;
			}

		VectorType expected(6);
		for (size_t col = 0; col < 2; ++col)
			for (size_t row = 0; row < 3; ++row)
				expected[row + col * 3] = Y(row, col);

		// Call matrixVector (implementation under test)
		bg.matrixVector(vout, vin);

		// Compare device result with host reference
		for (size_t i = 0; i < 6; ++i)
			CHECK(vout[i] == Catch::Approx(expected[i]).epsilon(1e-12));
	}

	if (useKokkos)
		Kokkos::finalize();
}

// TEST_CASE("BatchedGemmKokkos matrixVector", "[BatchedGemm]")
//{
//     runBatchedGemmTest<Dmrg::BatchedGemmKokkos>(true);
// }

TEST_CASE("BatchedGemmCpu matrixVector", "[BatchedGemm]")
{
	runBatchedGemmTest<Dmrg::BatchedGemmCpu>(true);
}

TEST_CASE("BatchedGemmPluginSc matrixVector", "[BatchedGemm]")
{
	runBatchedGemmTest<Dmrg::BatchedGemmPluginSc>(true);
}

TEST_CASE("BatchedGemm all implementations agree", "[BatchedGemm][integration]")
{
	Kokkos::initialize();
	FakeInitKron fk(1, 1);
	fk.lrs_.l.parts  = { 0, 2 };
	fk.lrs_.r.parts  = { 0, 3 };
	fk.patchLeft[0]  = 0;
	fk.patchRight[0] = 0;

	std::mt19937_64                        rng(123456);
	std::uniform_real_distribution<double> dist(-1.0, 1.0);

	PsimagLite::Matrix<double> A(2, 2), B(3, 3);
	for (int i = 0; i < 2; i++)
		for (int j = 0; j < 2; j++)
			A(i, j) = dist(rng);
	for (int i = 0; i < 3; i++)
		for (int j = 0; j < 3; j++)
			B(i, j) = dist(rng);
	fk.xc0.set(0, 0, FakeMatrixDenseOrSparse(A));
	fk.yc0.set(0, 0, FakeMatrixDenseOrSparse(B));

	using VectorType = PsimagLite::Vector<double>::Type;
	VectorType vin(6), v_k(6), v_c(6), v_p(6);
	for (size_t i = 0; i < 6; i++)
		vin[i] = dist(rng);

	// Expected host
	PsimagLite::Matrix<double> X(3, 2), BX(3, 2), Y(3, 2);
	for (size_t col = 0; col < 2; col++)
		for (size_t row = 0; row < 3; row++)
			X(row, col) = vin[row + col * 3];
	for (size_t i = 0; i < 3; i++)
		for (size_t j = 0; j < 2; j++) {
			double s = 0;
			for (size_t k = 0; k < 3; k++)
				s += B(i, k) * X(k, j);
			BX(i, j) = s;
		}
	for (size_t i = 0; i < 3; i++)
		for (size_t j = 0; j < 2; j++) {
			double s = 0;
			for (size_t k = 0; k < 2; k++)
				s += BX(i, k) * A(j, k);
			Y(i, j) = s;
		}
	VectorType expected(6);
	for (size_t col = 0; col < 2; col++)
		for (size_t row = 0; row < 3; row++)
			expected[row + col * 3] = Y(row, col);

	// Instantiate and run each implementation inside a scope so Kokkos objects are
	// destroyed before Kokkos::finalize()
	{
		// Dmrg::BatchedGemmKokkos<FakeInitKron> bgk(fk);
		Dmrg::BatchedGemmCpu<FakeInitKron>      bgc(fk);
		Dmrg::BatchedGemmPluginSc<FakeInitKron> bgp(fk);

		// bgk.matrixVector(v_k, vin);
		bgc.matrixVector(v_c, vin);
		bgp.matrixVector(v_p, vin);

		// Compare Kokkos and Plugin results immediately (they don't depend on
		// Kokkos::finalize)
		for (size_t i = 0; i < 6; i++) {
			// CHECK(v_k[i] == Catch::Approx(expected[i]).epsilon(1e-12));
			CHECK(v_p[i] == Catch::Approx(expected[i]).epsilon(1e-12));
			CHECK(v_c[i] == Catch::Approx(expected[i]).epsilon(1e-12));
		}
	}

	Kokkos::finalize();
}
