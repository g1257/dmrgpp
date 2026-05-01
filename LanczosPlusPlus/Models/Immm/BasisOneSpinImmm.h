
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

#ifndef BASIS_ONE_SPIN_IMMM_H
#define BASIS_ONE_SPIN_IMMM_H
#include "BitManip.h"
#include "LanczosPlusPlus/src/Engine/LabeledOperator.h"
#include "LanczosPlusPlus/src/Engine/LanczosGlobals.h"
#include "Matrix.h"
#include <cassert>

namespace LanczosPlusPlus {

class BasisOneSpinImmm {

public:

	typedef LanczosGlobals::WordType WordType;
	typedef LabeledOperator          LabeledOperatorType;

	BasisOneSpinImmm(const PsimagLite::Vector<SizeType>::Type& orbsPerSite, SizeType npart)
	    : orbsPerSite_(orbsPerSite)
	    , nsite_(orbsPerSite.size())
	    , npart_(npart)
	{
		LanczosGlobals::doCombinatorial(orbs() * nsite_ + 1);
		LanczosGlobals::doBitmask(nsite_ * 4 + 1);

		/* compute size of basis */
		if (npart == 0) {
			data_.resize(1);
			data_[0] = 0;
			return;
		}
		SizeType levels = 0;
		for (SizeType i = 0; i < orbsPerSite_.size(); i++)
			levels += orbsPerSite_[i];
		SizeType tmp = LanczosGlobals::combinatorial(levels, npart);
		data_.resize(tmp);

		levels = orbsPerSite_.size() * orbs();
		tmp    = LanczosGlobals::combinatorial(levels, npart);
		reordering_.resize(tmp);

		// compute basis:
		SizeType counter  = 0;
		SizeType counter2 = 0;
		for (SizeType na = 0; na <= npart; na++) {
			SizeType                           nb = npart - na;
			PsimagLite::Vector<WordType>::Type basisA, basisB;
			fillPartialBasis(basisA, na);
			fillPartialBasis(basisB, nb);
			collateBasis(counter, counter2, basisA, basisB);
		}
	}

	SizeType size() const { return data_.size(); }

	const WordType& operator[](SizeType i) const { return data_[i]; }

	void print(std::ostream& os, bool isBinary) const
	{
		if (isBinary) {
			LanczosGlobals::printBasisBinary(os, 1, data_);
		} else {
			LanczosGlobals::printBasisDecimal(os, 40, data_);
		}
	}

	SizeType perfectIndex(WordType ket) const
	{
		for (SizeType i = 0; i < data_.size(); i++)
			if (data_[i] == ket)
				return i;
		print(std::cout, true);
		assert(false);
		return 0;
	}

	SizeType getN(WordType ket, SizeType site, SizeType orb) const
	{
		WordType ketA = 0, ketB = 0;
		uncollateKet(ketA, ketB, ket);
		if (orb == 0) {
			WordType res = (ketA & LanczosGlobals::bitmask(site));
			return (res > 0) ? 1 : 0;
		}
		WordType res2 = ketB & LanczosGlobals::bitmask(site);
		return (res2 > 0) ? 1 : 0;
	}

	SizeType getN(SizeType i, SizeType orb) const
	{
		WordType ketA = 0, ketB = 0;
		uncollateKet(ketA, ketB, data_[i]);
		if (orb == 0)
			return PsimagLite::BitManip::count(ketA);
		return PsimagLite::BitManip::count(ketB);
	}

	SizeType getN(SizeType) const { throw std::runtime_error("getN\n"); }

	static const WordType& bitmask(SizeType i) { return LanczosGlobals::bitmask(i); }

	int doSign(SizeType i, SizeType site, SizeType orb) const
	{
		WordType ketA = 0, ketB = 0;
		uncollateKet(ketA, ketB, data_[i]);
		if (orb == 0) {
			return doSign(ketA, site);
		}

		SizeType c   = PsimagLite::BitManip::count(ketA);
		int      ret = (c & 1) ? LanczosGlobals::FERMION_SIGN : 1;
		return ret * doSign(ketB, site);
	}

	int doSign(WordType ket, SizeType i, SizeType orb1, SizeType j, SizeType orb2) const
	{
		if (i > j) {
			std::cerr << "FATAL: At doSign\n";
			std::cerr << "INFO: i=" << i << " j=" << j << std::endl;
			std::cerr << "AT: " << __FILE__ << " : " << __LINE__ << std::endl;
			throw std::runtime_error("BasisOneSpinImmm::doSign(...)\n");
		}
		SizeType x0 = (i + 1) * orbs(); // i+1 cannot be the last site, 'cause i<j
		SizeType x1 = j * orbs();

		SizeType sum = getNbyKet(ket, x0, x1);

		// at site i we need to be carefull
		x0 = i * orbs() + orb1;
		x1 = (i + 1) * orbs();
		sum += getNbyKet(ket, x0, x1);

		// same at site j
		x0 = j * orbs();
		x1 = j * orbs() + orb2;
		sum += getNbyKet(ket, x0, x1);

		return (sum & 1) ? LanczosGlobals::FERMION_SIGN : 1;
	}

	int doSignGf(WordType a, SizeType ind, SizeType orb) const
	{
		WordType ketA = 0, ketB = 0;

		uncollateKet(ketA, ketB, a);

		if (orb == 0)
			return doSignGf(ketA, ind);
		int s = (PsimagLite::BitManip::count(ketA) & 1) ? -1 : 1; // Parity of a

		return s * doSignGf(ketB, ind);
	}

	SizeType getNbyKet(SizeType ket) const
	{
		SizeType sum     = 0;
		WordType ketCopy = ket;
		while (ketCopy) {
			if (ketCopy & 1)
				sum++;
			ketCopy <<= 1;
		}
		return sum;
	}

	SizeType isThereAnElectronAt(SizeType ket, SizeType site, SizeType orb) const
	{
		SizeType x = site * orbs() + orb;
		return (ket & LanczosGlobals::bitmask(x)) ? 1 : 0;
	}

	SizeType electrons() const { return npart_; }

	bool getBra(WordType&                  bra,
	            const WordType&            myword,
	            const LabeledOperatorType& lOperator,
	            SizeType                   site,
	            SizeType                   orb) const
	{
		WordType ketA = 0, ketB = 0;
		uncollateKet(ketA, ketB, myword);
		WordType braA = ketA;
		WordType braB = ketB;

		if (orb == 0) {
			if (!getBra(braA, ketA, lOperator, site))
				return false;
		} else {
			if (!getBra(braB, ketB, lOperator, site))
				return false;
		}

		bra = getCollatedKet(braA, braB);
		return true;
	}

	int newPartCorCdagger(SizeType newPart1, const LabeledOperatorType& lOperator) const
	{
		int c = (lOperator.id() == LabeledOperatorType::Label::OPERATOR_CDAGGER) ? 1 : -1;
		newPart1 += c;

		if (SizeType(newPart1) > maxElectrons())
			return -1;

		return newPart1;
	}

private:

	bool getBra(WordType&                  bra,
	            const WordType&            ket,
	            const LabeledOperatorType& lOperator,
	            SizeType                   i) const
	{

		WordType si = (ket & LanczosGlobals::bitmask(i));
		if (lOperator.id() == LabeledOperatorType::Label::OPERATOR_C) {
			if (si > 0) {
				bra = (ket ^ LanczosGlobals::bitmask(i));
				return true;
			} else {
				return false; // cannot destroy, there's nothing
			}
		} else if (lOperator.id() == LabeledOperatorType::Label::OPERATOR_CDAGGER) {
			if (si == 0) {
				bra = (ket ^ LanczosGlobals::bitmask(i));
				return true;
			} else {
				return false; // cannot construct, there's already one
			}
		} else if (lOperator.id() == LabeledOperatorType::Label::OPERATOR_N) {
			if (si == 0)
				return false;
			bra = ket;
			return true;
		}

		PsimagLite::String str = lOperator.unknownOperator();
		throw std::runtime_error(str.c_str());
	}

	int doSignGf(WordType b, SizeType ind) const
	{
		if (ind == 0)
			return 1;

		// ind>0 from now on
		SizeType i    = 0;
		SizeType j    = ind;
		WordType mask = b;
		mask &= ((1 << (i + 1)) - 1) ^ ((1 << j) - 1);
		int s = (PsimagLite::BitManip::count(mask) & 1) ? -1
		                                                : 1; // Parity of up between i and j
		// Is there a down at i?
		if (LanczosGlobals::bitmask(i) & b)
			s = -s;
		return s;
	}

	SizeType maxElectrons() const
	{
		SizeType sum = 0;
		for (SizeType i = 0; i < orbsPerSite_.size(); i++)
			sum += orbsPerSite_[i];
		return sum;
	}

	SizeType orbs() const { return orbsPerSite_[0]; }

	void fillPartialBasis(PsimagLite::Vector<WordType>::Type& partialBasis, SizeType npart)
	{
		/* compute size of basis */
		SizeType hilbert = 1;
		int      n       = nsite_;
		SizeType m       = 1;
		for (; m <= npart; n--, m++)
			hilbert = hilbert * n / m;

		if (partialBasis.size() != hilbert) {
			partialBasis.clear();
			partialBasis.resize(hilbert);
		}

		if (npart == 0) {
			partialBasis[0] = 0;
			return;
		}
		/* define basis states */
		WordType ket = (1ul << npart) - 1;
		for (SizeType i = 0; i < hilbert; i++) {
			partialBasis[i] = ket;
			n = m = 0;
			for (; (ket & 3) != 1; n++, ket >>= 1) {
				m += ket & 1;
			}
			ket = ((ket + 1) << n) ^ ((1 << m) - 1);
		}
	}

	void collateBasis(SizeType&                                 counter,
	                  SizeType&                                 counter2,
	                  const PsimagLite::Vector<WordType>::Type& basisA,
	                  const PsimagLite::Vector<WordType>::Type& basisB)
	{
		for (SizeType i = 0; i < basisA.size(); i++) {
			for (SizeType j = 0; j < basisB.size(); j++) {
				reordering_[counter2++] = counter;
				if (isForbiddenSite(basisB[j]))
					continue;
				WordType ket = getCollatedKet(basisA[i], basisB[j]);
				assert(counter < data_.size());
				data_[counter++] = ket;
			}
		}
	}

	bool isForbiddenSite(const WordType& ket) const
	{
		for (SizeType i = 0; i < orbsPerSite_.size(); i++) {
			if (orbsPerSite_[i] > 1)
				continue;
			WordType mask = (1 << i);
			if (mask & ket)
				return true;
		}
		return false;
	}

	SizeType perfectIndexPartial(WordType state) const
	{
		SizeType n = 0;
		for (SizeType b = 0, c = 1; state > 0; b++, state >>= 1)
			if (state & 1)
				n += LanczosGlobals::combinatorial(b, c++);

		return n;
	}

	WordType getCollatedKet(WordType ketA, WordType ketB) const
	{
		WordType remA    = ketA;
		WordType remB    = ketB;
		SizeType counter = 0;
		WordType ket     = 0;

		while (remA || remB) {
			SizeType bitA = (remA & 1);
			SizeType bitB = (remB & 1);
			if (bitA)
				ket |= LanczosGlobals::bitmask(counter);
			if (bitB)
				ket |= LanczosGlobals::bitmask(counter + 1);
			counter += 2;
			if (remA)
				remA >>= 1;
			if (remB)
				remB >>= 1;
		}
		return ket;
	}

	void uncollateKet(WordType& ketA, WordType& ketB, WordType ket) const
	{
		SizeType counter = 0;
		ketA = ketB = 0;
		while (ket) {
			SizeType bitA = (ket & 1);
			SizeType bitB = (ket & 2);
			if (bitA)
				ketA |= LanczosGlobals::bitmask(counter);
			if (bitB)
				ketB |= LanczosGlobals::bitmask(counter);
			counter++;
			ket >>= 2;
		}
	}

	int doSign(WordType a, SizeType i) const
	{
		if (i == nsite_ - 1)
			return 1;

		a &= ((1 << (i + 1)) - 1) ^ ((1 << nsite_) - 1);
		// Parity of single occupied between i and nsite-1
		int s = (PsimagLite::BitManip::count(a) & 1) ? LanczosGlobals::FERMION_SIGN : 1;
		return s;
	}

	SizeType getNbyKet(SizeType ket, SizeType from, SizeType upto) const
	{
		SizeType sum     = 0;
		SizeType counter = from;
		while (counter < upto) {
			if (ket & LanczosGlobals::bitmask(counter))
				sum++;
			counter++;
		}
		return sum;
	}

	const PsimagLite::Vector<SizeType>::Type& orbsPerSite_;
	SizeType                                  nsite_;
	SizeType                                  npart_;
	PsimagLite::Vector<WordType>::Type        data_;
	PsimagLite::Vector<WordType>::Type        reordering_;

}; // class BasisOneSpinImmm

std::ostream& operator<<(std::ostream& os, const BasisOneSpinImmm& b);
} // namespace LanczosPlusPlus
#endif
