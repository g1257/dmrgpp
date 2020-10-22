#ifndef CINC_BASIS_H
#define CINC_BASIS_H
#include "BasisOneSpin.h"
#include "Vector.h"

namespace Dmft {

class BasisExactDiag {

public:

	typedef BasisOneSpin::WordType WordType;
	typedef BasisOneSpin::LabeledOperatorType LabeledOperatorType;

	BasisExactDiag(SizeType sites, SizeType nup,SizeType ndown)
	    : basis1_(sites, nup),
	      basis2_(sites, ndown)
	{}

	SizeType size() const { return basis1_.size()*basis2_.size(); }

	WordType operator()(SizeType i,SizeType spin) const
	{
		SizeType y = i/basis1_.size();
		SizeType x = i%basis1_.size();
		assert(x < basis1_.size());
		assert(y < basis2_.size());
		return (spin == 0) ? basis1_[x] : basis2_[y];
	}

	SizeType isThereAnElectronAt(WordType ket1,
	                             WordType ket2,
	                             SizeType site,
	                             SizeType spin,
	                             SizeType) const
	{
		return (spin==0) ? basis1_.isThereAnElectronAt(ket1,site)
		                 : basis2_.isThereAnElectronAt(ket2,site);
	}

	SizeType getN(WordType ket1,
	              WordType ket2,
	              SizeType site,
	              SizeType spin,
	              SizeType) const
	{
		return (spin == 0) ? basis1_.getN(ket1,site) : basis2_.getN(ket2,site);
	}

	static const WordType& bitmask(SizeType i)
	{
		return BasisOneSpin::bitmask(i);
	}

	SizeType perfectIndex(WordType ket1,WordType ket2) const
	{
		return basis1_.perfectIndex(ket1) +
		        basis2_.perfectIndex(ket2)*basis1_.size();
	}

	int doSignGf(WordType a,
	             WordType b,
	             SizeType ind,
	             SizeType sector,
	             SizeType) const
	{
		if (sector == 0) { // spin up
			if (ind==0) return 1;

			// ind>0 from now on
			SizeType i = 0;
			SizeType j = ind;
			WordType mask = a;
			mask &= ((1 << (i+1)) - 1) ^ ((1 << j) - 1);
			int s=(PsimagLite::BitManip::count(mask) & 1) ? -1 : 1;
			// Is there an up at i?
			if (bitmask(i) & a) s = -s;
			return s;
		}

		int s=(PsimagLite::BitManip::count(a) & 1) ? -1 : 1; // Parity of up
		if (ind==0) return s;

		// ind>0 from now on
		SizeType i = 0;
		SizeType j = ind;
		WordType mask = b;
		mask &= ((1 << (i+1)) - 1) ^ ((1 << j) - 1);
		s=(PsimagLite::BitManip::count(mask) & 1) ? -1 : 1;
		// Is there a down at i?
		if (bitmask(i) & b) s = -s;
		return s;
	}

	bool getBra(WordType& bra,
	            WordType ket1,
	            WordType ket2,
	            const LabeledOperatorType& lOperator,
	            SizeType site,
	            SizeType spin) const
	{
		return (spin == 0) ? basis1_.getBra(bra, ket1, lOperator, site) :
		                         basis2_.getBra(bra, ket2, lOperator, site);
	}


private:

	BasisOneSpin basis1_;
	BasisOneSpin basis2_;
};

}
#endif // CINC_BASIS_H
