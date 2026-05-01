
/*
// BEGIN LICENSE BLOCK
Copyright (c) 2009 , UT-Battelle, LLC
All rights reserved

[Lanczos++, Version 1.0.0]

*********************************************************
THE SOFTWARE IS SUPPLIED BY THE COPYRIGHT HOLDERS AND
CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED
WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A
PARTICULAR PURPOSE ARE DISCLAIMED.

Please see full open source license included in file LICENSE.
*********************************************************

*/

#ifndef BASIS_FEASBASED_SC_H
#define BASIS_FEASBASED_SC_H
#include "../../Engine/BasisBase.h"
#include "BasisOneSpinFeAs.h"

namespace LanczosPlusPlus {

template <typename GeometryType_> class BasisFeAsBasedSc : public BasisBase<GeometryType_> {

	typedef LanczosGlobals::PairIntType PairIntType;

	static const SizeType SPIN_UP = LanczosGlobals::SPIN_UP;

public:

	typedef GeometryType_                     GeometryType;
	typedef BasisBase<GeometryType>           BaseType;
	typedef typename BaseType::WordType       WordType;
	typedef typename BaseType::VectorWordType VectorWordType;
	typedef BasisOneSpinFeAs                  BasisType;
	typedef BasisType::LabeledOperatorType    LabeledOperatorType;

	BasisFeAsBasedSc(const GeometryType& geometry,
	                 SizeType            nup,
	                 SizeType            ndown,
	                 SizeType            orbitals)
	    : orbitals_(orbitals)
	    , nup_(nup)
	    , ndown_(ndown)
	    , basis1_(geometry.numberOfSites(), nup, orbitals)
	    , basis2_(geometry.numberOfSites(), ndown, orbitals)
	{ }

	PairIntType parts() const { return PairIntType(nup_, ndown_); }

	static const WordType& bitmask(SizeType i) { return BasisType::bitmask(i); }

	SizeType dofs() const { return 2 * orbitals_; }

	SizeType size() const { return basis1_.size() * basis2_.size(); }

	virtual SizeType hilbertOneSite(SizeType) const { return (1 << (orbitals_ * 2)); }

	WordType operator()(SizeType i, SizeType spin) const
	{
		SizeType y = i / basis1_.size();
		SizeType x = i % basis1_.size();
		return (spin == SPIN_UP) ? basis1_[x] : basis2_[y];
	}

	SizeType perfectIndex(const VectorWordType& kets) const
	{
		assert(kets.size() == 2);
		return perfectIndex(kets[0], kets[1]);
	}

	SizeType perfectIndex(WordType ket1, WordType ket2) const
	{
		return basis1_.perfectIndex(ket1) + basis2_.perfectIndex(ket2) * basis1_.size();
	}

	SizeType perfectIndex(WordType, SizeType, SizeType) const
	{
		throw PsimagLite::RuntimeError("perfectIndex\n");
	}

	SizeType getN(SizeType i, SizeType spin, SizeType orb) const
	{
		SizeType y = i / basis1_.size();
		SizeType x = i % basis1_.size();
		return (spin == SPIN_UP) ? basis1_.getN(x, orb) : basis2_.getN(y, orb);
	}

	SizeType
	getN(WordType ket1, WordType ket2, SizeType site, SizeType spin, SizeType orb) const
	{
		return (spin == SPIN_UP) ? basis1_.getN(ket1, site, orb)
		                         : basis2_.getN(ket2, site, orb);
	}

	PairIntType getBraIndex(WordType                   ket1,
	                        WordType                   ket2,
	                        const LabeledOperatorType& lOperator,
	                        SizeType                   site,
	                        SizeType                   spin,
	                        SizeType                   orb) const
	{
		if (lOperator.id() == LabeledOperatorType::Label::OPERATOR_C
		    || lOperator.id() == LabeledOperatorType::Label::OPERATOR_CDAGGER
		    || lOperator.id() == LabeledOperatorType::Label::OPERATOR_N)
			return PairIntType(
			    getBraIndexCorCdaggerOrN(ket1, ket2, lOperator, site, spin, orb), 1);

		if (lOperator.id() == LabeledOperatorType::Label::OPERATOR_SPLUS
		    || lOperator.id() == LabeledOperatorType::Label::OPERATOR_SMINUS)
			return PairIntType(
			    getBraIndexSplusOrSminus(ket1, ket2, lOperator, site, orb), 1);

		if (lOperator.id() == LabeledOperatorType::Label::OPERATOR_CDAGGER_A_UP_C_B_UP)
			return PairIntType(getBraIndexCdaggerC(ket1, ket2, site, 0, 0, 1), 1);

		PsimagLite::String str(__FILE__);
		str += " " + ttos(__LINE__) + "\n";
		str += PsimagLite::String("getBraIndex: unsupported operator ");
		str += lOperator.toString() + "\n";
		throw std::runtime_error(str.c_str());
	}

	int doSign(WordType ket1,
	           WordType ket2,
	           SizeType i,
	           SizeType orb1,
	           SizeType j,
	           SizeType orb2,
	           SizeType spin) const
	{
		if (i > j) {
			std::cerr << "FATAL: At doSign\n";
			std::cerr << "INFO: i=" << i << " j=" << j << std::endl;
			std::cerr << "AT: " << __FILE__ << " : " << __LINE__ << std::endl;
			throw std::runtime_error("FeBasedSc::doSign(...)\n");
		}
		if (spin == SPIN_UP) {
			return basis1_.doSign(ket1, i, orb1, j, orb2);
		}
		return basis2_.doSign(ket2, i, orb1, j, orb2);
	}

	int doSignGf(WordType a, WordType b, SizeType ind, SizeType spin, SizeType orb) const
	{
		if (spin == SPIN_UP)
			return basis1_.doSignGf(a, ind, orb);

		int s  = (PsimagLite::BitManip::count(a) & 1) ? -1 : 1; // Parity of up
		int s2 = basis2_.doSignGf(b, ind, orb);

		return s * s2;
	}

	int doSignSpinOrbit(WordType a,
	                    WordType b,
	                    SizeType ind,
	                    SizeType spin1,
	                    SizeType orb1,
	                    SizeType spin2,
	                    SizeType orb2) const
	{
		if (spin1 == spin2) {
			if (spin1 == 0)
				return basis1_.doSign(a, ind, orb1, orb2);
			return basis2_.doSign(b, ind, orb1, orb2);
		}

		int x = (spin1) ? -1 : 1;
		int s = (PsimagLite::BitManip::count(a) & 1) ? -1 : 1; // Parity of up
		if (spin1)
			return x * s * basis1_.doSignGf(a, ind, orb2)
			    * basis2_.doSignGf(b, ind, orb1);

		return x * s * basis1_.doSignGf(a, ind, orb1) * basis2_.doSignGf(b, ind, orb2);
	}

	int doSignSpSm(WordType a, WordType b, SizeType ind, SizeType spin, SizeType orb) const
	{
		if (spin == SPIN_UP) { // spin here means S^\dagger
			// FIXME: Count over a (up)
			return basis1_.doSignGf(a, ind, orb) * basis2_.doSignGf(b, ind, orb);
		}

		// FIXME: Count over a + 1
		return basis1_.doSignGf(a, ind, orb) * basis2_.doSignGf(b, ind, orb);
	}

	SizeType isThereAnElectronAt(WordType ket1,
	                             WordType ket2,
	                             SizeType site,
	                             SizeType spin,
	                             SizeType orb) const
	{
		if (spin == SPIN_UP)
			return basis1_.isThereAnElectronAt(ket1, site, orb);
		return basis2_.isThereAnElectronAt(ket2, site, orb);
	}

	bool hasNewParts(std::pair<SizeType, SizeType>&       newParts,
	                 const std::pair<SizeType, SizeType>& oldParts,
	                 const LabeledOperatorType&           lOperator,
	                 SizeType                             spin,
	                 SizeType                             orb) const
	{
		if (lOperator.id() == LabeledOperatorType::Label::OPERATOR_C
		    || lOperator.id() == LabeledOperatorType::Label::OPERATOR_CDAGGER)
			return hasNewPartsCorCdagger(newParts, oldParts, lOperator, spin);

		if (lOperator.id() == LabeledOperatorType::Label::OPERATOR_SPLUS
		    || lOperator.id() == LabeledOperatorType::Label::OPERATOR_SMINUS)
			return hasNewPartsSplusOrSminus(newParts, oldParts, lOperator);

		PsimagLite::String str(__FILE__);
		str += " " + ttos(__LINE__) + "\n";
		str += PsimagLite::String("hasNewParts: unsupported operator ");
		str += lOperator.toString() + "\n";
		throw std::runtime_error(str.c_str());
	}

	SizeType orbsPerSite(SizeType) const { throw PsimagLite::RuntimeError("orbsPerSite\n"); }

	SizeType orbs() const { throw PsimagLite::RuntimeError("orbs\n"); }

	void print(std::ostream& os, typename BaseType::PrintEnum binaryOrDecimal) const
	{
		bool isBinary = (binaryOrDecimal == BaseType::PRINT_BINARY);
		os << "\tUp sector\n";
		basis1_.print(os, isBinary);
		os << "\tDown sector\n";
		basis2_.print(os, isBinary);
	}

	virtual bool
	getBra(WordType&, WordType, WordType, const LabeledOperatorType&, SizeType, SizeType) const
	{
		throw PsimagLite::RuntimeError("BasisFeAs: getBra\n");
	}

private:

	int getBraIndexCorCdaggerOrN(WordType                   ket1,
	                             WordType                   ket2,
	                             const LabeledOperatorType& lOperator,
	                             SizeType                   site,
	                             SizeType                   spin,
	                             SizeType                   orb) const
	{

		WordType bra = 0;
		bool     b   = getBraCorCdaggerOrN(bra, ket1, ket2, lOperator, site, spin, orb);
		if (!b)
			return -1;
		return (spin == SPIN_UP) ? perfectIndex(bra, ket2) : perfectIndex(ket1, bra);
	}

	int getBraIndexSplusOrSminus(WordType                   ket1,
	                             WordType                   ket2,
	                             const LabeledOperatorType& lOperator,
	                             SizeType                   site,
	                             SizeType                   orb) const
	{

		WordType bra1 = 0;
		WordType bra2 = 0;
		bool     b    = getBraSplusOrSminus(bra1, bra2, ket1, ket2, lOperator, site, orb);
		if (!b)
			return -1;
		return perfectIndex(bra1, bra2);
	}

	bool hasNewPartsCorCdagger(std::pair<SizeType, SizeType>&       newParts,
	                           const std::pair<SizeType, SizeType>& oldParts,
	                           const LabeledOperatorType&           lOperator,
	                           SizeType                             spin) const
	{
		int newPart1 = oldParts.first;
		int newPart2 = oldParts.second;

		if (spin == SPIN_UP)
			newPart1 = basis1_.newPartCorCdagger(newPart1, lOperator);
		else
			newPart2 = basis2_.newPartCorCdagger(newPart2, lOperator);

		if (newPart1 < 0 || newPart2 < 0)
			return false;

		if (newPart1 == 0 && newPart2 == 0)
			return false;
		newParts.first  = SizeType(newPart1);
		newParts.second = SizeType(newPart2);
		return true;
	}

	bool hasNewPartsSplusOrSminus(std::pair<SizeType, SizeType>&       newParts,
	                              const std::pair<SizeType, SizeType>& oldParts,
	                              const LabeledOperatorType&           lOperator) const
	{
		int c1 = (lOperator.id() == LabeledOperatorType::Label::OPERATOR_SPLUS) ? 1 : -1;
		int c2 = (lOperator.id() == LabeledOperatorType::Label::OPERATOR_SPLUS) ? -1 : 1;

		int newPart1 = basis1_.hasNewPartsSplusOrSminus(oldParts.first, c1);
		int newPart2 = basis2_.hasNewPartsSplusOrSminus(oldParts.second, c2);

		if (newPart1 < 0 || newPart2 < 0)
			return false;

		if (newPart1 == 0 && newPart2 == 0)
			return false;
		newParts.first  = SizeType(newPart1);
		newParts.second = SizeType(newPart2);
		return true;
	}

	bool getBraCorCdaggerOrN(WordType&                  bra,
	                         const WordType&            ket1,
	                         const WordType&            ket2,
	                         const LabeledOperatorType& lOperator,
	                         SizeType                   site,
	                         SizeType                   spin,
	                         SizeType                   orb) const
	{
		return (spin == SPIN_UP) ? basis1_.getBra(bra, ket1, lOperator, site, orb)
		                         : basis2_.getBra(bra, ket2, lOperator, site, orb);
	}

	bool getBraSplusOrSminus(WordType&                  bra1,
	                         WordType&                  bra2,
	                         const WordType&            ket1,
	                         const WordType&            ket2,
	                         const LabeledOperatorType& lOperator,
	                         SizeType                   site,
	                         SizeType                   orb) const
	{
		const LabeledOperatorType::Label what1
		    = (lOperator.id() == LabeledOperatorType::Label::OPERATOR_SPLUS)
		    ? LabeledOperatorType::Label::OPERATOR_CDAGGER
		    : LabeledOperatorType::Label::OPERATOR_C;

		const LabeledOperatorType::Label what2
		    = (lOperator.id() == LabeledOperatorType::Label::OPERATOR_SPLUS)
		    ? LabeledOperatorType::Label::OPERATOR_C
		    : LabeledOperatorType::Label::OPERATOR_CDAGGER;

		bool b1 = basis1_.getBra(bra1, ket1, LabeledOperatorType(what1), site, orb);

		bool b2 = basis2_.getBra(bra2, ket2, LabeledOperatorType(what2), site, orb);

		return (b1 & b2);
	}

	int getBraIndexCdaggerC(WordType ket1,
	                        WordType ket2,
	                        SizeType site,
	                        SizeType spin,
	                        SizeType orb1,
	                        SizeType orb2) const
	{
		assert(spin == 0 || spin == 1);
		const BasisType& mybasis = (spin == 0) ? basis1_ : basis2_;
		WordType         ket     = (spin == 0) ? ket1 : ket2;
		WordType         bra     = 0;
		bool             b1      = mybasis.getBra(
                    bra, ket, LabeledOperatorType::Label::OPERATOR_CDAGGER, site, orb1);

		if (!b1)
			return -1;

		WordType newbra = 0;
		bool     b2     = mybasis.getBra(
                    newbra, bra, LabeledOperatorType::Label::OPERATOR_C, site, orb2);
		if (!b2)
			return -1;

		return (spin == 0) ? perfectIndex(newbra, ket2) : perfectIndex(ket1, newbra);
	}

	SizeType  orbitals_;
	SizeType  nup_;
	SizeType  ndown_;
	BasisType basis1_;
	BasisType basis2_;
}; // class BasisFeAsBasedSc
} // namespace LanczosPlusPlus
#endif
