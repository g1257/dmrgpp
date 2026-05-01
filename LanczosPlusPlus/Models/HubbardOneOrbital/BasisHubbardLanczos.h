/*
 */

#ifndef BASISHUBBARDLANCZOS_H
#define BASISHUBBARDLANCZOS_H
#include "../../Engine/BasisBase.h"
#include "BasisOneSpin.h"

namespace LanczosPlusPlus {

template <typename GeometryType> class BasisHubbardLanczos : public BasisBase<GeometryType> {

	enum
	{
		SPIN_UP   = LanczosGlobals::SPIN_UP,
		SPIN_DOWN = LanczosGlobals::SPIN_DOWN
	};

public:

	typedef LanczosGlobals::PairIntType            PairIntType;
	typedef BasisOneSpin                           BasisType;
	typedef BasisBase<GeometryType>                BaseType;
	typedef typename BaseType::WordType            WordType;
	typedef typename BaseType::VectorWordType      VectorWordType;
	typedef typename BaseType::LabeledOperatorType LabeledOperatorType;

	BasisHubbardLanczos(const GeometryType& geometry, SizeType nup, SizeType ndown)
	    : nup_(nup)
	    , ndown_(ndown)
	    , basis1_(geometry.numberOfSites(), nup)
	    , basis2_(geometry.numberOfSites(), ndown)
	{ }

	PairIntType parts() const { return PairIntType(nup_, ndown_); }

	static const WordType& bitmask(SizeType i) { return BasisType::bitmask(i); }

	SizeType size() const { return basis1_.size() * basis2_.size(); }

	//! Spin up and spin down
	SizeType dofs() const { return 2; }

	virtual SizeType hilbertOneSite(SizeType) const { return 4; }

	SizeType perfectIndex(const VectorWordType& kets) const
	{
		assert(kets.size() == 2);
		return perfectIndex(kets[0], kets[1]);
	}

	SizeType perfectIndex(WordType ket1, WordType ket2) const
	{
		return basis1_.perfectIndex(ket1) + basis2_.perfectIndex(ket2) * basis1_.size();
	}

	virtual SizeType perfectIndex(WordType, SizeType, SizeType) const
	{
		throw PsimagLite::RuntimeError("perfectIndex\n");
	}

	SizeType electrons(SizeType what) const
	{
		return (what == SPIN_UP) ? basis1_.electrons() : basis2_.electrons();
	}

	WordType operator()(SizeType i, SizeType spin) const
	{
		SizeType y = i / basis1_.size();
		SizeType x = i % basis1_.size();
		assert(x < basis1_.size());
		assert(y < basis2_.size());
		return (spin == SPIN_UP) ? basis1_[x] : basis2_[y];
	}

	SizeType isThereAnElectronAt(WordType ket1,
	                             WordType ket2,
	                             SizeType site,
	                             SizeType spin,
	                             SizeType) const
	{
		if (spin == SPIN_UP)
			return basis1_.isThereAnElectronAt(ket1, site);
		return basis2_.isThereAnElectronAt(ket2, site);
	}

	SizeType getN(WordType ket1, WordType ket2, SizeType site, SizeType spin, SizeType) const
	{
		return (spin == SPIN_UP) ? basis1_.getN(ket1, site) : basis2_.getN(ket2, site);
	}

	int doSignGf(WordType a, WordType b, SizeType ind, SizeType sector, SizeType) const
	{
		if (sector == SPIN_UP) {
			if (ind == 0)
				return 1;

			// ind>0 from now on
			SizeType i    = 0;
			SizeType j    = ind;
			WordType mask = a;
			mask &= ((1 << (i + 1)) - 1) ^ ((1 << j) - 1);
			int s = (PsimagLite::BitManip::count(mask) & 1)
			    ? LanczosGlobals::FERMION_SIGN
			    : 1;
			// Is there an up at i?
			if (BasisType::bitmask(i) & a)
				s = -s;
			return s;
		}
		int s = (PsimagLite::BitManip::count(a) & 1) ? LanczosGlobals::FERMION_SIGN
		                                             : 1; // Parity of up
		if (ind == 0)
			return s;

		// ind>0 from now on
		SizeType i    = 0;
		SizeType j    = ind;
		WordType mask = b;
		mask &= ((1 << (i + 1)) - 1) ^ ((1 << j) - 1);
		s = (PsimagLite::BitManip::count(mask) & 1) ? LanczosGlobals::FERMION_SIGN : 1;
		// Is there a down at i?
		if (BasisType::bitmask(i) & b)
			s = -s;
		return s;
	}

	int doSign(WordType ket1,
	           WordType ket2,
	           SizeType i,
	           SizeType,
	           SizeType j,
	           SizeType,
	           SizeType spin) const
	{
		assert(i <= j);
		return (spin == SPIN_UP) ? basis1_.doSign(ket1, i, j) : basis2_.doSign(ket2, i, j);
	}

	int doSignSpSm(WordType a, WordType b, SizeType ind, SizeType spin, SizeType) const
	{
		if (spin == SPIN_UP) { // spin here means S^\dagger
			// FIXME: Count over a (up)
			return LanczosGlobals::doSign(a, ind) * LanczosGlobals::doSign(b, ind);
		}

		// FIXME: Count over a + 1
		return LanczosGlobals::doSign(a, ind) * LanczosGlobals::doSign(b, ind);
	}

	PairIntType getBraIndex(WordType                   ket1,
	                        WordType                   ket2,
	                        const LabeledOperatorType& lOperator,
	                        SizeType                   site,
	                        SizeType                   spin,
	                        SizeType) const
	{
		if (lOperator.id() == LabeledOperatorType::Label::OPERATOR_SPLUS
		    || lOperator.id() == LabeledOperatorType::Label::OPERATOR_SMINUS)
			return getBraIndexSplusSminus(ket1, ket2, lOperator, site);

		if (lOperator.id() == LabeledOperatorType::Label::OPERATOR_SZ)
			return getBraIndexSz(ket1, ket2, site);

		WordType bra = 0;
		bool     b   = getBra(bra, ket1, ket2, lOperator, site, spin);
		if (!b)
			return PairIntType(-1, 1);
		int tmp = (spin == SPIN_UP) ? perfectIndex(bra, ket2) : perfectIndex(ket1, bra);
		return PairIntType(tmp, 1);
	}

	SizeType orbsPerSite(SizeType) const { return 1; }

	SizeType orbs() const { return 1; }

	void print(std::ostream& os, typename BaseType::PrintEnum binaryOrDecimal) const
	{
		bool isBinary = (binaryOrDecimal == BaseType::PRINT_BINARY);
		os << "\tUp sector\n";
		basis1_.print(os, isBinary);
		os << "\tDown sector\n";
		basis2_.print(os, isBinary);
	}

private:

	bool getBra(WordType&                  bra,
	            WordType                   ket1,
	            WordType                   ket2,
	            const LabeledOperatorType& lOperator,
	            SizeType                   site,
	            SizeType                   spin) const
	{
		return (spin == SPIN_UP) ? basis1_.getBra(bra, ket1, lOperator, site)
		                         : basis2_.getBra(bra, ket2, lOperator, site);
	}

	PairIntType getBraIndexSz(WordType ket1, WordType ket2, SizeType site) const
	{
		LabeledOperatorType opN(LabeledOperatorType::Label::OPERATOR_N);
		WordType            bra = 0;
		bool                b1  = basis1_.getBra(bra, ket1, opN, site);
		bool                b2  = basis2_.getBra(bra, ket2, opN, site);
		if (!b1 && !b2)
			return PairIntType(-1, 1);
		if (b1 && b2)
			return PairIntType(-1, 1);
		int      tmp   = (b1) ? 1 : -1;
		SizeType index = perfectIndex(ket1, ket2);
		return PairIntType(index, tmp);
	}

	PairIntType getBraIndexSplusSminus(WordType                   ket1,
	                                   WordType                   ket2,
	                                   const LabeledOperatorType& lOperator,
	                                   SizeType                   site) const
	{
		SizeType spin = (lOperator.id() == LabeledOperatorType::Label::OPERATOR_SPLUS)
		    ? SPIN_UP
		    : SPIN_DOWN;

		LabeledOperatorType opC(LabeledOperatorType::Label::OPERATOR_C);
		WordType            brar1 = 0;
		bool b = getBra(brar1, ket1, ket2, opC.transposeConjugate(), site, spin);
		if (!b)
			return PairIntType(-1, 1);

		WordType brar2 = 0;
		b              = getBra(brar2, ket1, ket2, opC, site, 1 - spin);
		if (!b)
			return PairIntType(-1, 1);

		int tmp
		    = (spin == SPIN_UP) ? perfectIndex(brar1, brar2) : perfectIndex(brar2, brar1);

		return PairIntType(tmp, 1);
	}

	SizeType  nup_;
	SizeType  ndown_;
	BasisType basis1_;
	BasisType basis2_;

}; // class BasisHubbardLanczos

} // namespace LanczosPlusPlus
#endif
