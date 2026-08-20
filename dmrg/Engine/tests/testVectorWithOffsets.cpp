#include "Qn.h"
#include "VectorWithOffsets.h"
#include <catch2/catch_approx.hpp>
#include <catch2/catch_template_test_macros.hpp>
#include <catch2/catch_test_macros.hpp>
#include <complex>
#include <vector>

namespace {

using VectorSizeType = std::vector<SizeType>;

template <typename T> void checkValue(const T& actual, double real, double imag = 0.0)
{
	CHECK(std::real(actual) == Catch::Approx(real));
	CHECK(std::imag(actual) == Catch::Approx(imag));
}

Dmrg::Qn makeQn(bool odd)
{
	return Dmrg::Qn(odd, VectorSizeType(), Dmrg::Qn::PairSizeType(0, 0), 0);
}

class FakeBasis {

public:

	FakeBasis()
	    : partitions_({ 0, 2, 5, 6 })
	    , qns_({ makeQn(false), makeQn(true), makeQn(false) })
	{ }

	SizeType        size() const { return partitions_.back(); }
	SizeType        partition() const { return partitions_.size(); }
	SizeType        partition(SizeType sector) const { return partitions_[sector]; }
	const Dmrg::Qn& pseudoQn(SizeType sector) const { return qns_[sector]; }
	const Dmrg::Qn& qnEx(SizeType sector) const { return qns_[sector]; }

private:

	VectorSizeType        partitions_;
	std::vector<Dmrg::Qn> qns_;
};

} // namespace

TEMPLATE_TEST_CASE("VectorWithOffsets default state and identity",
                   "[VectorWithOffsets]",
                   double,
                   std::complex<double>)
{
	using VectorWithOffsetsType = Dmrg::VectorWithOffsets<TestType, Dmrg::Qn>;
	VectorWithOffsetsType v;
	CHECK(v.size() == 0);
	CHECK(v.sectors() == 0);
	CHECK(VectorWithOffsetsType::name() == "vectorwithoffsets");
	CHECK(norm(v) == 0.0);
}

TEMPLATE_TEST_CASE("VectorWithOffsets constructs populated sectors",
                   "[VectorWithOffsets]",
                   double,
                   std::complex<double>)
{
	using VectorWithOffsetsType = Dmrg::VectorWithOffsets<TestType, Dmrg::Qn>;
	const FakeBasis basis;

	SECTION("weights for all basis sectors")
	{
		VectorSizeType        weights({ 2, 0, 1 });
		VectorWithOffsetsType v(weights, basis);

		CHECK(v.size() == 6);
		CHECK(v.sectors() == 2);
		CHECK(v.sector(0) == 0);
		CHECK(v.sector(1) == 2);
		CHECK(v.qn(0) == basis.pseudoQn(0));
		CHECK(v.qn(1) == basis.pseudoQn(2));
		CHECK(v.offset(0) == 0);
		CHECK(v.offset(1) == 2);
		CHECK(v.offset(2) == 5);
		CHECK(v.offset(3) == 6);
		CHECK(v.effectiveSize(0) == 2);
		CHECK(v.effectiveSize(1) == 0);
		CHECK(v.effectiveSize(2) == 1);
		CHECK(v.index2Sector(0) == 0);
		CHECK(v.index2Sector(2) == -1);
		CHECK(v.index2Sector(5) == 2);
	}

	SECTION("compacted weights and explicit sectors")
	{
		VectorSizeType        weights({ 3, 1 });
		VectorSizeType        sectors({ 1, 2 });
		VectorWithOffsetsType v(weights, sectors, basis);
		CHECK(v.sectors() == 2);
		CHECK(v.sector(0) == 1);
		CHECK(v.sector(1) == 2);
		CHECK(v.effectiveSize(0) == 0);
		CHECK(v.effectiveSize(1) == 3);
		CHECK(v.effectiveSize(2) == 1);
	}
}

TEMPLATE_TEST_CASE("VectorWithOffsets sets, extracts, and sparsifies data",
                   "[VectorWithOffsets]",
                   double,
                   std::complex<double>)
{
	using VectorType            = std::vector<TestType>;
	using VectorWithOffsetsType = Dmrg::VectorWithOffsets<TestType, Dmrg::Qn>;
	const FakeBasis       basis;
	VectorType            data({ 3.0, 4.0, 5.0 });
	VectorWithOffsetsType v;
	v.set(data, 1, basis);

	CHECK(data.empty());
	CHECK(v.size() == 6);
	CHECK(v.sectors() == 1);
	CHECK(v.sector(0) == 1);
	CHECK(v.qn(0) == basis.pseudoQn(1));
	CHECK(v.index2Sector(1) == -1);
	CHECK(v.index2Sector(2) == 1);

	VectorType extracted;
	v.extract(extracted, 1);
	CHECK(extracted == VectorType({ 3.0, 4.0, 5.0 }));

	VectorType replacement({ 6.0, 7.0, 8.0 });
	v.setDataInSector(replacement, 1);
	checkValue(v.fastAccess(1, 0), 6.0);
	v.fastAccess(1, 1) = 9.0;
	checkValue(v.slowAccess(3), 9.0);
	v.slowAccess(4) = 10.0;

	const VectorWithOffsetsType& constV = v;
	checkValue(constV.slowAccess(0), 0.0);
	checkValue(constV.slowAccess(4), 10.0);

	std::vector<TestType> sparse;
	v.toSparse(sparse);
	CHECK(sparse == std::vector<TestType>({ 0.0, 0.0, 6.0, 9.0, 10.0, 0.0 }));

	v.clear();
	CHECK(v.size() == 0);
	CHECK(v.sectors() == 0);
}

TEMPLATE_TEST_CASE("VectorWithOffsets converts a full vector",
                   "[VectorWithOffsets]",
                   double,
                   std::complex<double>)
{
	using VectorType            = std::vector<TestType>;
	using VectorWithOffsetsType = Dmrg::VectorWithOffsets<TestType, Dmrg::Qn>;
	const FakeBasis       basis;
	VectorWithOffsetsType v;
	v.fromFull(VectorType({ 1.0, 2.0, 0.0, 0.0, 0.0, 3.0 }), basis);

	CHECK(v.size() == 6);
	CHECK(v.sectors() == 2);
	CHECK(v.sector(0) == 0);
	CHECK(v.sector(1) == 2);
	CHECK(v.fastAccess(0, 0) == 1.0);
	CHECK(v.fastAccess(0, 1) == 2.0);
	CHECK(v.fastAccess(2, 0) == 3.0);
	CHECK(v.index2Sector(3) == -1);
}

TEMPLATE_TEST_CASE("VectorWithOffsets populates and collapses sectors",
                   "[VectorWithOffsets]",
                   double,
                   std::complex<double>)
{
	using VectorWithOffsetsType = Dmrg::VectorWithOffsets<TestType, Dmrg::Qn>;
	const FakeBasis       basis;
	VectorWithOffsetsType v;
	v.populateSectors(basis);
	CHECK(v.sectors() == 3);
	CHECK(v.effectiveSize(0) == 2);
	CHECK(v.effectiveSize(1) == 3);
	CHECK(v.effectiveSize(2) == 1);

	v.fastAccess(1, 1) = 2.0;
	v.collapseSectors();
	CHECK(v.sectors() == 1);
	CHECK(v.sector(0) == 1);
	CHECK(v.index2Sector(0) == -1);
	CHECK(v.index2Sector(2) == 1);

	VectorWithOffsetsType copy;
	copy.populateFromQns(v, basis);
	CHECK(copy.sectors() == 1);
	CHECK(copy.sector(0) == 1);
	CHECK(copy.effectiveSize(1) == 3);
	CHECK(copy.fastAccess(1, 0) == 0.0);
}

TEMPLATE_TEST_CASE("VectorWithOffsets public arithmetic",
                   "[VectorWithOffsets]",
                   double,
                   std::complex<double>)
{
	using VectorType            = std::vector<TestType>;
	using VectorWithOffsetsType = Dmrg::VectorWithOffsets<TestType, Dmrg::Qn>;
	const FakeBasis       basis;
	VectorType            data({ 3.0, 4.0, 0.0 });
	VectorWithOffsetsType v;
	v.set(data, 1, basis);

	CHECK(norm(v) == Catch::Approx(5.0));
	normalize(v);
	CHECK(norm(v) == Catch::Approx(1.0));

	v *= 2.0;
	checkValue(v.fastAccess(1, 0), 1.2);
	VectorWithOffsetsType scaled = 0.5 * v;
	checkValue(scaled.fastAccess(1, 1), 0.8);
	checkValue(scaled * scaled, 1.0);

	VectorWithOffsetsType sum;
	sum += scaled;
	sum += scaled;
	checkValue(sum.fastAccess(1, 0), 1.2);
	checkValue(sum.fastAccess(1, 1), 1.6);
}
