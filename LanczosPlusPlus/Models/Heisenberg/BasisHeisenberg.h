/*
 */

#ifndef LANCZOS_BASIS_HEISENBERG_H
#define LANCZOS_BASIS_HEISENBERG_H

#include "BasisBase.h"
#include "BitManip.h"
#include "LanczosPlusPlus/src/Engine/LanczosGlobals.h"

namespace LanczosPlusPlus {

template <typename GeometryType> class BasisHeisenberg : public BasisBase<GeometryType> {

public:

	typedef LanczosGlobals::PairIntType            PairIntType;
	typedef BasisBase<GeometryType>                BaseType;
	typedef typename BaseType::WordType            WordType;
	typedef typename BaseType::VectorWordType      VectorWordType;
	typedef typename BaseType::LabeledOperatorType LabeledOperatorType;

	BasisHeisenberg(const GeometryType& geometry, SizeType twiceS, SizeType szPlusConst)
	    : geometry_(geometry)
	    , twiceS_(twiceS)
	    , szPlusConst_(szPlusConst)
	    , bits_(0)
	{
		SizeType sites = geometry_.numberOfSites();
		LanczosGlobals::doBitmask(sites);
		WordType searchTotal = 1;
		assert(twiceS > 0);
		bits_ = 1 + static_cast<SizeType>(logBase2(twiceS + 1));
		if (twiceS & 1)
			bits_--;
		searchTotal <<= (bits_ * sites);

		WordType mask = getMask();

		for (WordType lui = 0; lui < searchTotal; ++lui) {
			int tmp = mOf(lui, mask);
			if (tmp < 0 || static_cast<SizeType>(tmp) != szPlusConst)
				continue;
			data_.push_back(lui);
		}
	}

	PairIntType parts() const { return PairIntType(twiceS_, szPlusConst_); }

	static const WordType& bitmask(SizeType i) { return LanczosGlobals::bitmask(i); }

	SizeType size() const { return data_.size(); }

	SizeType dofs() const { return twiceS_ + 1; }

	virtual SizeType hilbertOneSite(SizeType) const { return 1 + twiceS_; }

	SizeType perfectIndex(const VectorWordType&) const
	{
		throw PsimagLite::RuntimeError("BasisHeisenberg::perfectIndex kets\n");
	}

	SizeType perfectIndex(WordType ket, WordType) const
	{
		for (SizeType i = 0; i < data_.size(); ++i) {
			if (ket == data_[i])
				return i;
		}

		throw PsimagLite::RuntimeError("perfectIndex: no index found\n");
	}

	WordType operator()(SizeType i, SizeType) const { return data_[i]; }

	SizeType isThereAnElectronAt(WordType, WordType, SizeType, SizeType, SizeType) const
	{
		throw PsimagLite::RuntimeError("BasisHeisenberg::isThereAnElectronAt\n");
	}

	SizeType getN(WordType ket1, WordType, SizeType site, SizeType, SizeType) const
	{
		WordType mask = getMask();
		ket1 >>= (bits_ * site);
		return (ket1 & mask);
	}

	int doSignGf(WordType, WordType, SizeType, SizeType, SizeType) const
	{
		throw PsimagLite::RuntimeError("BasisHeisenberg::doSignGf\n");
	}

	int doSign(WordType, WordType, SizeType, SizeType, SizeType, SizeType, SizeType) const
	{
		throw PsimagLite::RuntimeError("BasisHeisenberg::doSign\n");
	}

	PairIntType getBraIndex(WordType                   ket1,
	                        WordType                   ket2,
	                        const LabeledOperatorType& lOperator,
	                        SizeType                   site,
	                        SizeType                   spin,
	                        SizeType                   orb) const
	{
		if (twiceS_ != 1)
			throw PsimagLite::RuntimeError("BasisHeisenberg::getBraIndex_ \n");

		if (lOperator.id() == LabeledOperatorType::Label::OPERATOR_SPLUS
		    || lOperator.id() == LabeledOperatorType::Label::OPERATOR_SMINUS) {
			return getBraIndexSplusSminus(ket1, ket2, lOperator, site, spin, orb);
		}

		return getBraIndex_(ket1, ket2, lOperator, site, spin, orb);
	}

	SizeType orbsPerSite(SizeType) const { return 1; }

	SizeType orbs() const { return 1; }

	SizeType perfectIndex(WordType, SizeType, SizeType) const
	{
		throw PsimagLite::RuntimeError("BasisHeisenberg::perfectIndex\n");
	}

	void print(std::ostream& os, typename BaseType::PrintEnum binaryOrDecimal) const
	{
		SizeType hilbert = 1;
		hilbert <<= geometry_.numberOfSites();
		if (binaryOrDecimal == BaseType::PRINT_BINARY) {
			LanczosGlobals::printBasisBinary(os, hilbert, data_);
		} else {
			LanczosGlobals::printBasisDecimal(os, 40, data_);
		}
	}

	SizeType szPlusConst() const { return szPlusConst_; }

	template <typename GeometryType2>
	friend std::ostream& operator<<(std::ostream&                         os,
	                                const BasisHeisenberg<GeometryType2>& basis);

private:

	bool getBra(WordType&                  bra,
	            WordType                   ket,
	            WordType                   site1,
	            const LabeledOperatorType& val1,
	            SizeType                   site2,
	            SizeType                   val2) const
	{
		bra            = ket;
		WordType mask1 = getMask();
		WordType mask2 = mask1;
		mask1 <<= (site1 * bits_);
		bra &= (~mask1);

		mask2 <<= (site2 * bits_);
		bra &= (~mask2);

		mask1 = val1.toUint();
		mask1 <<= (site1 * bits_);
		bra |= mask1;

		mask2 = val2;
		mask2 <<= (site2 * bits_);
		bra |= mask2;

		return true;
	}

	WordType getMask() const
	{
		SizeType mask = 1;
		for (SizeType i = 0; i < bits_; ++i)
			mask |= LanczosGlobals::bitmask(i);
		return mask;
	}

	int mOf(WordType lui, WordType mask) const
	{
		SizeType m       = 0;
		bool     allowed = true;
		while (lui != 0) {
			WordType tmp = (lui & mask);
			if (!mOfIsAllowed(tmp)) {
				allowed = false;
				break;
			}

			m += tmp;
			lui >>= bits_;
		}

		return (allowed) ? m : -1;
	}

	bool mOfIsAllowed(WordType val) const
	{
		if (twiceS_ & 1)
			return true;

		return (val <= twiceS_);
	}

	PairIntType getBraIndexSplusSminus(WordType                   ket1,
	                                   WordType                   ket2,
	                                   const LabeledOperatorType& lOperator,
	                                   SizeType                   site,
	                                   SizeType,
	                                   SizeType) const
	{
		assert(lOperator.id() == LabeledOperatorType::Label::OPERATOR_SPLUS
		       || lOperator.id() == LabeledOperatorType::Label::OPERATOR_SMINUS);

		WordType mask1 = getMask();
		mask1 <<= (site * bits_);

		if (ket1 & mask1) {
			if (lOperator.id() == LabeledOperatorType::Label::OPERATOR_SPLUS)
				return PairIntType(-1, 0);
		} else {
			if (lOperator.id() == LabeledOperatorType::Label::OPERATOR_SMINUS)
				return PairIntType(-1, 0);
		}

		WordType bra = ket1 ^ mask1;

		return PairIntType(perfectIndex(bra, ket2), 1);
	}

	PairIntType getBraIndex_(WordType                   ket1,
	                         WordType                   ket2,
	                         const LabeledOperatorType& lOperator,
	                         SizeType                   site,
	                         SizeType                   spin,
	                         SizeType                   orb) const
	{
		if (twiceS_ != 1)
			err("getBraIndex_: S=1/2 supported only\n");

		if (lOperator.id() == LabeledOperatorType::Label::OPERATOR_N) {

			SizeType nup = getN(ket1, ket2, site, 0, orb);
			assert(nup < 2);
			return PairIntType(perfectIndex(ket1, ket2),
			                   (spin == LanczosGlobals::SPIN_UP) ? nup : 1 - nup);
		} else if (lOperator.id() == LabeledOperatorType::Label::OPERATOR_SZ) {
			int nup = getN(ket1, ket2, site, 0, orb);
			return PairIntType(perfectIndex(ket1, ket2), 1 - 2 * nup);
		}

		throw PsimagLite::RuntimeError(lOperator.unknownOperator());
	}

	SizeType logBase2(SizeType x) const
	{
		int ret = 0;
		while (x >>= 1)
			++ret;
		return ret;
	}

	const GeometryType& geometry_;
	SizeType            twiceS_;
	SizeType            szPlusConst_;
	SizeType            bits_;
	VectorWordType      data_;
}; // class BasisHeisenberg

template <typename GeometryType>
std::ostream& operator<<(std::ostream& os, const BasisHeisenberg<GeometryType>& basis)
{
	for (SizeType i = 0; i < basis.data_.size(); i++)
		os << i << " " << basis.data_[i] << "\n";
	return os;
}

} // namespace LanczosPlusPlus
#endif
