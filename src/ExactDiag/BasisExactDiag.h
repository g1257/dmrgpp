#ifndef CINC_BASIS_H
#define CINC_BASIS_H
#include "BasisOneSpin.h"
#include "Vector.h"

namespace Dmft {

class BasisExactDiag {

public:

	typedef BasisOneSpin::WordType WordType;

	BasisExactDiag(SizeType sites, SizeType nup,SizeType ndown)
	    : nup_(nup),
	      ndown_(ndown),
	      basis1_(sites, nup),
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

private:

	SizeType nup_;
	SizeType ndown_;
	BasisOneSpin basis1_;
	BasisOneSpin basis2_;
};

}
#endif // CINC_BASIS_H
