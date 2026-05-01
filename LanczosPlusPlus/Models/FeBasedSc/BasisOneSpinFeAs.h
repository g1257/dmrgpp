
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

#ifndef BASIS_ONE_SPIN_FE_AS_H
#define BASIS_ONE_SPIN_FE_AS_H
#include "BitManip.h"
#include "LabeledOperator.h"
#include "LanczosGlobals.h"
#include "Matrix.h"
#include "Partitions.h"

namespace LanczosPlusPlus {

class BasisOneSpinFeAs {

	typedef Partitions PartitionsType;

public:

	typedef LanczosGlobals::WordType WordType;
	typedef LabeledOperator          LabeledOperatorType;

	static PsimagLite::Matrix<SizeType> comb_;

	BasisOneSpinFeAs(SizeType nsite, SizeType npart, SizeType orbitals)
	    : orbitals_(orbitals)
	    , nsite_(nsite)
	    , npart_(npart)
	{
		LanczosGlobals::doCombinatorial(orbitals_ * nsite_ + 1);
		LanczosGlobals::doBitmask(nsite_ * orbitals_ * orbitals_ + 1);

		/* compute size of basis */
		if (npart == 0) {
			data_.resize(1);
			size_    = 1;
			data_[0] = 0;
			return;
		}
		size_ = 0;
		PartitionsType partitions(npart, orbitals_);
		for (SizeType i = 0; i < partitions.size(); i++) {
			const PsimagLite::Vector<SizeType>::Type& na  = partitions(i);
			SizeType                                  tmp = 1;
			for (SizeType j = 0; j < na.size(); j++)
				tmp *= LanczosGlobals::combinatorial(nsite_, na[j]);
			size_ += tmp;
		}
		data_.resize(size_);

		// compute basis:
		SizeType counter = 0;

		for (SizeType i = 0; i < partitions.size(); i++) {
			const PsimagLite::Vector<SizeType>::Type& na = partitions(i);
			PsimagLite::Vector<PsimagLite::Vector<WordType>::Type>::Type basisA(
			    orbitals_);
			for (SizeType orb = 0; orb < orbitals_; orb++) {
				fillPartialBasis(basisA[orb], na[orb]);
			}
			collateBasis(counter, basisA);
		}
	}

	SizeType size() const { return size_; }

	const WordType& operator[](SizeType i) const { return data_[i]; }

	SizeType perfectIndex(WordType ket) const
	{
		for (SizeType i = 0; i < data_.size(); i++)
			if (data_[i] == ket)
				return i;
		throw std::runtime_error("perfectindex\n");
	}

	SizeType getN(WordType ket, SizeType site, SizeType orb) const
	{
		PsimagLite::Vector<WordType>::Type kets(orbitals_, 0);
		uncollateKet(kets, ket);

		WordType res = (kets[orb] & LanczosGlobals::bitmask(site));
		return (res > 0) ? 1 : 0;
	}

	SizeType getN(SizeType i, SizeType orb) const
	{
		PsimagLite::Vector<WordType>::Type kets(orbitals_, 0);
		uncollateKet(kets, data_[i]);
		return PsimagLite::BitManip::count(kets[orb]);
	}

	SizeType getN(SizeType i) const
	{
		SizeType c = 0;
		for (SizeType orb = 0; orb < orbitals_; orb++)
			c += getN(i, orb);
		return c;
	}

	bool getBra(WordType&                  bra,
	            const WordType&            myword,
	            const LabeledOperatorType& lOperator,
	            SizeType                   site,
	            SizeType                   orb) const
	{
		PsimagLite::Vector<WordType>::Type kets(orbitals_, 0);
		uncollateKet(kets, myword);

		WordType braA = kets[orb];
		if (!getBraCorCdagger(braA, kets[orb], lOperator, site))
			return false;

		kets[orb] = braA;
		bra       = getCollatedKet(kets);
		return true;
	}

	static const WordType& bitmask(SizeType i) { return LanczosGlobals::bitmask(i); }

	int doSign(WordType ket, SizeType i, SizeType orb1, SizeType j, SizeType orb2) const
	{
		if (i > j) {
			std::cerr << "FATAL: At doSign\n";
			std::cerr << "INFO: i=" << i << " j=" << j << std::endl;
			std::cerr << "AT: " << __FILE__ << " : " << __LINE__ << std::endl;
			throw std::runtime_error("FeBasedSc::doSign(...)\n");
		}

		if (i == j)
			return doSign(ket, i, orb1, orb2);

		SizeType x0 = (i + 1) * orbitals_; // i+1 cannot be the last site, 'cause i<j
		SizeType x1 = j * orbitals_;

		SizeType sum = getNbyKet(ket, x0, x1);

		// at site i we need to be carefull
		x0 = i * orbitals_ + orb1;
		x1 = (i + 1) * orbitals_;
		sum += getNbyKet(ket, x0, x1);

		// same at site j
		x0 = j * orbitals_;
		x1 = j * orbitals_ + orb2;
		sum += getNbyKet(ket, x0, x1);

		return (sum & 1) ? LanczosGlobals::FERMION_SIGN : 1;
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
		SizeType x = site * orbitals_ + orb;
		return (ket & LanczosGlobals::bitmask(x)) ? 1 : 0;
	}

	SizeType electrons() const { return npart_; }

	int newPartCorCdagger(SizeType newPart1, const LabeledOperatorType& lOperator) const
	{
		int c = (lOperator.id() == LabeledOperatorType::Label::OPERATOR_C) ? -1 : 1;
		newPart1 += c;

		if (SizeType(newPart1) > orbitals_ * nsite_)
			return -1;

		return newPart1;
	}

	int hasNewPartsSplusOrSminus(int newPart1, int c) const
	{
		newPart1 += c;

		if (newPart1 < 0)
			return -1;

		if (SizeType(newPart1) > orbitals_ * nsite_)
			return -1;

		return newPart1;
	}

	int doSignGf(WordType a, SizeType ind, SizeType orb) const
	{
		SizeType x0  = 0;
		SizeType x1  = ind * orbitals_;
		SizeType sum = getNbyKet(a, x0, x1);

		// at site ind we need to be carefull
		x0 = ind * orbitals_;
		x1 = ind * orbitals_ + orb;
		sum += getNbyKet(a, x0, x1);

		return (sum & 1) ? LanczosGlobals::FERMION_SIGN : 1;
	}

	void print(std::ostream& os, bool isBinary) const
	{
		SizeType hilbert = 1;
		hilbert <<= (orbitals_ * nsite_);
		if (isBinary) {
			LanczosGlobals::printBasisBinary(os, hilbert, data_);
		} else {
			LanczosGlobals::printBasisDecimal(os, 40, data_);
		}
	}

	int doSign(WordType ket, SizeType i, SizeType orb1, SizeType orb2) const
	{
		if (orb1 > orb2)
			return -doSign(ket, i, orb2, orb1);

		SizeType x0  = i * orbitals_ + orb1;
		SizeType x1  = i * orbitals_ + orb2;
		SizeType sum = getNbyKet(ket, x0, x1);
		return (sum & 1) ? LanczosGlobals::FERMION_SIGN : 1;
	}

private:

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

	void
	collateBasis(SizeType&                                                           counter,
	             const PsimagLite::Vector<PsimagLite::Vector<WordType>::Type>::Type& basisA)
	{
		SizeType total = 1;
		for (SizeType orb = 0; orb < orbitals_; orb++) {
			total *= basisA[orb].size();
		}

		for (SizeType i = 0; i < total; i++) {
			PsimagLite::Vector<WordType>::Type kets(orbitals_);
			getKets(kets, basisA, i);
			WordType ket     = getCollatedKet(kets);
			data_[counter++] = ket;
		}
	}

	void getKets(PsimagLite::Vector<WordType>::Type&                                 kets,
	             const PsimagLite::Vector<PsimagLite::Vector<WordType>::Type>::Type& basisA,
	             SizeType                                                            ind) const
	{
		// ind = i0 + i1 * size0 + i2 * size0 * size1 + ...
		SizeType tmp = ind;

		SizeType sizes = 1;
		for (SizeType orb = 0; orb < orbitals_ - 1; orb++)
			sizes *= basisA[orb].size();

		for (SizeType orb = 1; orb < orbitals_; orb++) {
			SizeType ix           = SizeType(tmp / sizes);
			tmp                   = ind % sizes;
			kets[orbitals_ - orb] = basisA[orbitals_ - orb][ix];
			sizes /= basisA[orbitals_ - orb - 1].size();
		}
		kets[0] = basisA[0][tmp];
	}

	void doCombinatorial()
	{
		/* look-up table for binomial coefficients */
		comb_.resize(orbitals_ * nsite_ + 1, orbitals_ * nsite_ + 1, 0);

		for (SizeType n = 0; n < comb_.n_row(); n++) {
			SizeType m   = 0;
			int      j   = n;
			SizeType i   = 1;
			SizeType cnm = 1;
			for (; m <= n / 2; m++, cnm = cnm * j / i, i++, j--)
				comb_(n, m) = comb_(n, n - m) = cnm;
		}
	}

	SizeType perfectIndexPartial(WordType state) const
	{
		SizeType n = 0;
		for (SizeType b = 0, c = 1; state > 0; b++, state >>= 1)
			if (state & 1)
				n += comb_(b, c++);

		return n;
	}

	WordType getCollatedKet(const PsimagLite::Vector<WordType>::Type& kets) const
	{
		SizeType counter = 0;
		WordType ket     = 0;

		PsimagLite::Vector<WordType>::Type remA = kets;

		while (orAll(remA)) {
			for (SizeType orb = 0; orb < kets.size(); orb++) {
				SizeType bitA = (remA[orb] & 1);
				if (bitA)
					ket |= LanczosGlobals::bitmask(counter);
				counter++;
				if (remA[orb])
					remA[orb] >>= 1;
			}
		}

		return ket;
	}

	SizeType orAll(const PsimagLite::Vector<WordType>::Type& kets) const
	{
		SizeType b = 0;
		for (SizeType orb = 0; orb < kets.size(); orb++)
			b |= kets[orb];
		return b;
	}

	//! kets must be all zero here
	void uncollateKet(PsimagLite::Vector<WordType>::Type& kets, WordType ket) const
	{
		SizeType counter = 0;

		while (ket) {
			for (SizeType orb = 0; orb < kets.size(); orb++) {
				SizeType mask = (1 << orb);
				SizeType bitA = (ket & mask);

				if (bitA)
					kets[orb] |= LanczosGlobals::bitmask(counter);
			}
			counter++;
			ket >>= orbitals_;
		}
	}

	bool getBraCorCdagger(WordType&                  bra,
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

	SizeType                           orbitals_;
	SizeType                           nsite_;
	SizeType                           size_;
	SizeType                           npart_;
	PsimagLite::Vector<WordType>::Type data_;
}; // class BasisOneSpinFeAs

std::ostream& operator<<(std::ostream& os, const BasisOneSpinFeAs& b);

} // namespace LanczosPlusPlus
#endif
