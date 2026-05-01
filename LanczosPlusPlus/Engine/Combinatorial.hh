#ifndef COMBINATORIAL_HH
#define COMBINATORIAL_HH
#include "Matrix.h"

namespace LanczosPlusPlus {

class Combinatorial {
public:

	explicit Combinatorial(SizeType m)
	    : comb_(m, m)
	{
		doCombinatorial();
	}

	Combinatorial() { }

	void resize(SizeType m)
	{
		comb_.resize(m, m, 0);
		doCombinatorial();
	}

	SizeType size() const
	{
		assert(comb_.rows() == comb_.cols());
		return comb_.rows();
	}

	// PsimagLite::Matrix does assert for out-of-range
	const SizeType& operator()(SizeType i, SizeType j) const { return comb_(i, j); }

private:

	void doCombinatorial()
	{
		/* look-up table for binomial coefficients */
		for (SizeType n = 0; n < comb_.rows(); ++n) {
			SizeType m   = 0;
			int      j   = n;
			SizeType i   = 1;
			SizeType cnm = 1;
			for (; m <= n / 2; m++, cnm = cnm * j / i, i++, j--)
				comb_(n, m) = comb_(n, n - m) = cnm;
		}
	}

	PsimagLite::Matrix<SizeType> comb_;
};
}
#endif // COMBINATORIAL_HH
