/*
 */

#ifndef LANCZOSPP_BASIS_FERMIONSPINLESS_HH
#define LANCZOSPP_BASIS_FERMIONSPINLESS_HH
#include "../../Engine/BasisBase.h"
#include "../HubbardOneOrbital/BasisOneSpin.h"

namespace LanczosPlusPlus {

template <typename GeometryType> class BasisFermionSpinless : public BasisBase<GeometryType> {

public:

	typedef LanczosGlobals::PairIntType            PairIntType;
	typedef BasisOneSpin                           BasisType;
	typedef BasisBase<GeometryType>                BaseType;
	typedef typename BaseType::WordType            WordType;
	typedef typename BaseType::VectorWordType      VectorWordType;
	typedef typename BaseType::LabeledOperatorType LabeledOperatorType;

	BasisFermionSpinless(const GeometryType& geometry, SizeType ne)
	    : ne_(ne)
	    , basis_(geometry.numberOfSites(), ne)
	{ }

	PairIntType parts() const override { return PairIntType(ne_, 0); }

	static const WordType& bitmask(SizeType i) { return BasisType::bitmask(i); }

	SizeType size() const override { return basis_.size(); }

	SizeType dofs() const override { return 1; }

	SizeType hilbertOneSite(SizeType) const override { return 2; }

	SizeType perfectIndex(const VectorWordType& kets) const override
	{
		assert(kets.size() == 1);
		return perfectIndex(kets[0], 0);
	}

	SizeType perfectIndex(WordType ket1, WordType) const override
	{
		return basis_.perfectIndex(ket1);
	}

	SizeType perfectIndex(WordType, SizeType, SizeType) const override
	{
		throw PsimagLite::RuntimeError("perfectIndex\n");
	}

	WordType operator()(SizeType i, SizeType) const override
	{
		assert(i < basis_.size());
		return basis_[i];
	}

	SizeType isThereAnElectronAt(WordType ket,
	                             WordType,
	                             SizeType site,
	                             SizeType,
	                             SizeType) const override
	{
		return basis_.isThereAnElectronAt(ket, site);
	}

	SizeType getN(WordType ket, WordType, SizeType site, SizeType, SizeType) const override
	{
		return basis_.getN(ket, site);
	}

	int doSignGf(WordType a, WordType b, SizeType ind, SizeType, SizeType) const override
	{
		if (ind == 0)
			return 1;

		// ind>0 from now on
		SizeType i    = 0;
		SizeType j    = ind;
		WordType mask = a;
		mask &= ((1 << (i + 1)) - 1) ^ ((1 << j) - 1);
		int s = (PsimagLite::BitManip::count(mask) & 1) ? -1 : 1;
		// Is there an electron at i?
		if (BasisType::bitmask(i) & a)
			s = -s;
		return s;
	}

	int doSign(WordType ket1, WordType, SizeType i, SizeType, SizeType j, SizeType, SizeType)
	    const override
	{
		assert(i <= j);
		return basis_.doSign(ket1, i, j);
	}

	int doSignSpSm(WordType a, WordType b, SizeType ind, SizeType spin, SizeType) const override
	{
		throw PsimagLite::RuntimeError("doSignSpSm\n");
	}

	PairIntType getBraIndex(WordType ket1,
	                        WordType,
	                        const LabeledOperatorType& lOperator,
	                        SizeType                   site,
	                        SizeType,
	                        SizeType) const override
	{
		WordType bra = 0;
		bool     b   = getBra(bra, ket1, 0, lOperator, site, 0);
		if (!b)
			return PairIntType(-1, 1);
		int tmp = perfectIndex(bra, 0);
		return PairIntType(tmp, 1);
	}

	SizeType orbsPerSite(SizeType) const override { return 1; }

	SizeType orbs() const override { return 1; }

	void print(std::ostream& os, typename BaseType::PrintEnum binaryOrDecimal) const override
	{
		bool isBinary = (binaryOrDecimal == BaseType::PRINT_BINARY);
		basis_.print(os, isBinary);
	}

	SizeType electrons() const { return ne_; }

private:

	bool getBra(WordType& bra,
	            WordType  ket1,
	            WordType,
	            const LabeledOperatorType& lOperator,
	            SizeType                   site,
	            SizeType) const override
	{
		return basis_.getBra(bra, ket1, lOperator, site);
	}

	SizeType  ne_;
	BasisType basis_;

}; // class BasisFermionSpinless

} // namespace LanczosPlusPlus
#endif
