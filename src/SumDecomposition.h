#ifndef PSI_SUM_DECOMPOSITION_H
#define PSI_SUM_DECOMPOSITION_H

#include <iostream>
#include <cassert>
#include "Vector.h"
#include <cstdlib>
#include <numeric>

namespace PsimagLite {

class SumDecomposition {

	typedef Vector<SizeType>::Type VectorSizeType;

public:

	enum SelEnum {SEL_ALL, SEL_SIZE, SEL_INDEX};

	SumDecomposition(SizeType total,
	                 SizeType sum,
	                 SelEnum sel = SEL_ALL,
	                 int selection = 0)
	    : sum_(sum),sel_(sel),selection_(selection),size_(0)
	{
		VectorSizeType data(total,0);
		reentrant1_(data,sum,0);
	}

	SizeType size() const { return size_; }

	SizeType storedSize() const { return data_.size(); }

	const VectorSizeType& operator()(SizeType i) const
	{
		if (i < data_.size()) return data_[i];
		String msg("SumDecomposition::operator():");
		throw RuntimeError(msg + " requested index is too big\n");
	}

private:

	void reentrant1_(VectorSizeType& data,SizeType n, SizeType x)
	{
		SizeType sum = 0;
		if (x >= data.size()) {
			if (std::accumulate(data.begin(),data.end(),sum) != sum_) return;
			bool b1 = (sel_ == SEL_ALL);
			bool b2 = (sel_ == SEL_INDEX && selection_ == size_);
			if (b1 || b2) data_.push_back(data);
			size_++;
			return;
		}

		for (SizeType i = 0; i <= n; ++i) {
			data[x] = i;
			reentrant1_(data,n-i,x+1);
		}
	}

	SizeType sum_;
	SelEnum sel_;
	int selection_;
	int size_;
	Vector<VectorSizeType>::Type data_;

};

std::ostream& operator<<(std::ostream& os, const SumDecomposition& sd)
{
	os<<sd.storedSize()<<"\n";
	for (SizeType i = 0; i < sd.storedSize(); ++i)
		os<<sd(i)<<"\n";

	return os;
}

} // namespace PsimagLite

#endif // PSI_SUM_DECOMPOSITION_H

