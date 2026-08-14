#include "Qn.h"
#include "VectorWithOffset.h"
#include <catch2/catch_approx.hpp>
#include <catch2/catch_template_test_macros.hpp>
#include <catch2/catch_test_macros.hpp>
#include <complex>
#include <map>
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
	    : partitions_({ 0, 2, 5 })
	    , qns_({ makeQn(false), makeQn(true) })
	{ }

	SizeType size() const { return partitions_.back(); }

	SizeType partition() const { return partitions_.size(); }

	SizeType partition(SizeType sector) const { return partitions_[sector]; }

	const Dmrg::Qn& pseudoQn(SizeType sector) const { return qns_[sector]; }

	const Dmrg::Qn& qnEx(SizeType sector) const { return qns_[sector]; }

private:

	VectorSizeType        partitions_;
	std::vector<Dmrg::Qn> qns_;
};

template <typename T> class FakeIo {

public:

	using VectorType           = std::vector<T>;
	using VectorWithOffsetType = Dmrg::VectorWithOffset<T, Dmrg::Qn>;
	using PairQnType           = typename VectorWithOffsetType::PairQnType;

	void createGroup(const PsimagLite::String& label) { group_ = label; }

	void write(SizeType value, const PsimagLite::String& label) { sizes_[label] = value; }

	void write(const PairQnType& value, const PsimagLite::String& label)
	{
		pairs_.insert_or_assign(label, value);
	}

	void write(const VectorType& value, const PsimagLite::String& label)
	{
		vectors_[label] = value;
	}

	void read(SizeType& value, const PsimagLite::String& label) { value = sizes_.at(label); }

	void read(PairQnType& value, const PsimagLite::String& label) { value = pairs_.at(label); }

	void read(VectorType& value, const PsimagLite::String& label)
	{
		value = vectors_.at(label);
	}

	const PsimagLite::String& group() const { return group_; }

private:

	PsimagLite::String                       group_;
	std::map<PsimagLite::String, SizeType>   sizes_;
	std::map<PsimagLite::String, PairQnType> pairs_;
	std::map<PsimagLite::String, VectorType> vectors_;
};

} // namespace

TEMPLATE_TEST_CASE("VectorWithOffset default state and metadata",
                   "[VectorWithOffset]",
                   double,
                   std::complex<double>)
{
	using VectorWithOffsetType = Dmrg::VectorWithOffset<TestType, Dmrg::Qn>;
	VectorWithOffsetType v;

	CHECK(v.size() == 0);
	CHECK(v.sectors() == 0);
	CHECK(v.effectiveSize() == 0);
	CHECK(v.offset() == 0);
	CHECK(v.index2Sector(0) == -1);
	CHECK(VectorWithOffsetType::name() == "vectorwithoffset");
	CHECK(VectorWithOffsetType::zero_ == TestType(0.0));
}

TEMPLATE_TEST_CASE("VectorWithOffset constructors expose one sector",
                   "[VectorWithOffset]",
                   double,
                   std::complex<double>)
{
	using VectorWithOffsetType = Dmrg::VectorWithOffset<TestType, Dmrg::Qn>;
	const FakeBasis basis;

	SECTION("weight and sector constructor")
	{
		VectorWithOffsetType v(3, 1, basis);
		CHECK(v.size() == 5);
		CHECK(v.sectors() == 1);
		CHECK(v.sector(0) == 1);
		CHECK(v.qn(0) == basis.pseudoQn(1));
		CHECK(v.offset(0) == 2);
		CHECK(v.effectiveSize(0) == 3);
	}

	SECTION("compacted weights constructor")
	{
		VectorSizeType       weights(1, 2);
		VectorSizeType       sectors(1, 0);
		VectorWithOffsetType v(weights, sectors, basis);
		CHECK(v.size() == 5);
		CHECK(v.sector(0) == 0);
		CHECK(v.offset() == 0);
		CHECK(v.effectiveSize() == 2);
	}

	SECTION("multiple sectors are rejected")
	{
		VectorSizeType weights(2, 1);
		VectorSizeType sectors(2, 0);
		CHECK_THROWS_AS(VectorWithOffsetType(weights, sectors, basis),
		                PsimagLite::RuntimeError);
	}
}

TEMPLATE_TEST_CASE("VectorWithOffset sets, extracts, and clears data",
                   "[VectorWithOffset]",
                   double,
                   std::complex<double>)
{
	using VectorType           = std::vector<TestType>;
	using VectorWithOffsetType = Dmrg::VectorWithOffset<TestType, Dmrg::Qn>;
	const FakeBasis      basis;
	VectorType           data({ 1.0, 2.0, 3.0 });
	VectorWithOffsetType v;
	v.set(data, 1, basis);

	CHECK(data.empty());
	CHECK(v.size() == 5);
	CHECK(v.sector(0) == 1);
	CHECK(v.qn(0) == basis.pseudoQn(1));
	CHECK(v.offset() == 2);

	VectorType extracted;
	v.extract(extracted);
	CHECK(extracted == VectorType({ 1.0, 2.0, 3.0 }));

	VectorType replacement({ 4.0, 5.0, 6.0 });
	v.setDataInSector(replacement, 0);
	v.extract(extracted, 0);
	CHECK(extracted == replacement);

	std::vector<TestType> sparse;
	v.toSparse(sparse);
	CHECK(sparse == std::vector<TestType>({ 0.0, 0.0, 4.0, 5.0, 6.0 }));

	v.clear();
	CHECK(v.size() == 0);
	CHECK(v.sectors() == 0);
	CHECK(v.effectiveSize() == 0);
	CHECK(v.offset() == 0);
}

TEMPLATE_TEST_CASE("VectorWithOffset converts from a full vector",
                   "[VectorWithOffset]",
                   double,
                   std::complex<double>)
{
	using VectorType           = std::vector<TestType>;
	using VectorWithOffsetType = Dmrg::VectorWithOffset<TestType, Dmrg::Qn>;
	const FakeBasis      basis;
	VectorWithOffsetType v;
	v.fromFull(VectorType({ 0.0, 0.0, 7.0, 8.0, 9.0 }), basis);

	CHECK(v.size() == 5);
	CHECK(v.sector(0) == 1);
	CHECK(v.qn(0) == basis.pseudoQn(1));
	CHECK(v.offset() == 2);
	CHECK(v.effectiveSize() == 3);
	CHECK(v.fastAccess(0, 0) == 7.0);
	CHECK(v.slowAccess(4) == 9.0);
	CHECK(v.index2Sector(1) == -1);
	CHECK(v.index2Sector(2) == 0);
	CHECK(v.index2Sector(4) == 0);
	CHECK(v.index2Sector(5) == -1);
}

TEMPLATE_TEST_CASE("VectorWithOffset public arithmetic",
                   "[VectorWithOffset]",
                   double,
                   std::complex<double>)
{
	using VectorType           = std::vector<TestType>;
	using VectorWithOffsetType = Dmrg::VectorWithOffset<TestType, Dmrg::Qn>;
	const FakeBasis      basis;
	VectorType           first({ 3.0, 4.0, 0.0 });
	VectorWithOffsetType v;
	v.set(first, 1, basis);

	CHECK(norm(v) == Catch::Approx(5.0));
	normalize(v);
	checkValue(v.fastAccess(0, 0), 0.6);
	checkValue(v.fastAccess(0, 1), 0.8);

	v.fastAccess(0, 2) = 1.0;
	v.slowAccess(2)    = 1.0;
	v *= 2.0;
	checkValue(v.fastAccess(0, 0), 2.0);

	VectorWithOffsetType scaled = 0.5 * v;
	checkValue(scaled.fastAccess(0, 0), 1.0);
	checkValue(scaled * scaled, 2.64);

	VectorWithOffsetType sum;
	sum += scaled;
	sum += scaled;
	checkValue(sum.fastAccess(0, 0), 2.0);
	checkValue(sum.fastAccess(0, 1), 1.6);
	checkValue(sum.fastAccess(0, 2), 2.0);

	VectorWithOffsetType otherSector(2, 0, basis);
	CHECK(scaled * otherSector == TestType(0.0));
	CHECK_THROWS_AS(sum += otherSector, PsimagLite::RuntimeError);
}

TEMPLATE_TEST_CASE("VectorWithOffset populates matching quantum-number sector",
                   "[VectorWithOffset]",
                   double,
                   std::complex<double>)
{
	using VectorType           = std::vector<TestType>;
	using VectorWithOffsetType = Dmrg::VectorWithOffset<TestType, Dmrg::Qn>;
	const FakeBasis      basis;
	VectorType           data({ 1.0, 2.0, 3.0 });
	VectorWithOffsetType source;
	source.set(data, 1, basis);

	VectorWithOffsetType destination;
	destination.populateFromQns(source, basis);
	CHECK(destination.size() == 5);
	CHECK(destination.sector(0) == 1);
	CHECK(destination.qn(0) == basis.pseudoQn(1));
	CHECK(destination.offset() == 2);
	CHECK(destination.effectiveSize() == 3);
	CHECK(destination.fastAccess(0, 0) == 0.0);

	destination.collapseSectors();
	CHECK(destination.sectors() == 1);
	CHECK_THROWS_AS(destination.populateSectors(basis), PsimagLite::RuntimeError);
}

TEMPLATE_TEST_CASE("VectorWithOffset writes and reads its public representation",
                   "[VectorWithOffset]",
                   double,
                   std::complex<double>)
{
	using VectorType           = std::vector<TestType>;
	using VectorWithOffsetType = Dmrg::VectorWithOffset<TestType, Dmrg::Qn>;
	const FakeBasis      basis;
	VectorType           data({ 2.0, 4.0, 6.0 });
	VectorWithOffsetType original;
	original.set(data, 1, basis);

	FakeIo<TestType> io;
	original.write(io, "vector");
	CHECK(io.group() == "vector");

	VectorWithOffsetType restored;
	restored.read(io, "vector");
	CHECK(restored.size() == original.size());
	CHECK(restored.offset() == original.offset());
	CHECK(restored.sector(0) == original.sector(0));
	CHECK(restored.qn(0) == original.qn(0));
	checkValue(restored * original, 56.0);

	VectorWithOffsetType loaded;
	loaded.loadOneSector(io, "vector", 7);
	checkValue(loaded * original, 56.0);
}
