#ifndef CINC_BASIS_ONE_SPIN_H
#define CINC_BASIS_ONE_SPIN_H
#include "Vector.h"
#include "Matrix.h"

namespace Dmft {

class BasisOneSpin {

public:

	typedef long unsigned int WordType;

	BasisOneSpin(SizeType nsite, SizeType npart)
	    : npart_(npart)
	{
		if (nsite_>0 && nsite!=nsite_)
			err("BasisOneSpin: All basis must have same number of sites\n");
		nsite_ = nsite;
		doCombinatorial();
		doBitmask(nsite_);

		/* compute size of basis */
		SizeType hilbert=1;
		int n=nsite;
		SizeType m=1;
		for (;m<=npart;n--,m++)
			hilbert=hilbert*n/m;

		if (data_.size()!=hilbert) {
			data_.clear();
			data_.resize(hilbert);
		}

		if (npart==0) {
			data_[0]=0;
			size_ = 1;
			return;
		}

		/* define basis states */
		WordType ket = (1ul<<npart)-1;
		for (SizeType i=0;i<hilbert;i++) {
			data_[i] = ket;
			n=m=0;
			for (;(ket&3)!=1;n++,ket>>=1) {
				m += ket&1;
			}
			ket = ((ket+1)<<n) ^ ((1<<m)-1);
		}
		size_ = hilbert;
	}

	SizeType size() const { return size_; }

	const WordType& operator[](SizeType i) const
	{
		return data_[i];
	}

	static SizeType isThereAnElectronAt(WordType ket,SizeType site)
	{
		return (ket & bitmask(site)) ? 1 : 0;
	}

	static SizeType getN(WordType ket,SizeType site)
	{
		return isThereAnElectronAt(ket,site);
	}

	static const WordType& bitmask(SizeType i)
	{
		assert(i < bitmask_.size());
		return bitmask_[i];
	}

	SizeType perfectIndex(WordType state) const
	{
		SizeType n=0;
		for (SizeType b=0,c=1;state>0;b++,state>>=1)
			if (state&1) n += comb_(b,c++);

		assert(n<data_.size());
		return n;
	}

private:

	static void doCombinatorial()
	{
		/* look-up table for binomial coefficients */
		comb_.resize(2*nsite_ + 2, 2*nsite_ + 2, 0);

		for (SizeType n=0;n<comb_.n_row();n++) {
			SizeType m = 0;
			int j = n;
			SizeType i = 1;
			SizeType cnm  = 1;
			for (;m<=n/2;m++,cnm=cnm*j/i,i++,j--)
				comb_(n,m) = comb_(n,n-m) = cnm;
		}
	}

	static void doBitmask(SizeType total)
	{
		if (total == bitmask_.size())
			return;

		bitmask_.resize(total);
		bitmask_[0]=1ul;
		for (SizeType i=1; i<bitmask_.size(); i++)
			bitmask_[i] = bitmask_[i-1]<<1;
	}

	static SizeType nsite_;
	static PsimagLite::Matrix<SizeType> comb_;
	static PsimagLite::Vector<WordType>::Type bitmask_;
	SizeType size_;
	SizeType npart_;
	PsimagLite::Vector<WordType>::Type data_;
};

}
#endif
