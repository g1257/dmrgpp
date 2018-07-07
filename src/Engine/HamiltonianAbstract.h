#ifndef HAMILTONIANABSTRACT_H
#define HAMILTONIANABSTRACT_H

#include "Vector.h"

namespace Dmrg {

class HamiltonianAbstract {

	typedef PsimagLite::Vector<SizeType>::Type VectorSizeType;
	typedef PsimagLite::Vector<VectorSizeType>::Type VectorVectorSizeType;

public:

	HamiltonianAbstract(const VectorSizeType& block)
	    : block_(block), data_(block.size()*block.size())
	{
		VectorSizeType v(2, 0);
		SizeType n = block.size();
		SizeType counter = 0;
		for (SizeType i = 0; i < n; ++i) {
			for (SizeType j = 0; j < n; ++j) {
				v[0] = block[i];
				v[1] = block[j];
				data_[counter++] = v;
			}
		}

		assert(counter == n*n);
	}

	SizeType items() const { return data_.size(); }

	VectorSizeType item(SizeType ind) const
	{
		assert(ind < data_.size());
		return data_[ind];
	}

	const VectorSizeType& block() const { return block_; }

private:

	const VectorSizeType& block_;
	VectorVectorSizeType data_;
};
}
#endif // HAMILTONIANABSTRACT_H
