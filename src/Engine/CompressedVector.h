#ifndef COMPRESSEDVECTOR_H
#define COMPRESSEDVECTOR_H
#include "Vector.h"
#include <fstream>

namespace Dmrg {

template<typename T>
class CompressedVector {

public:

	typedef T value_type;

	CompressedVector()
	    : size_(0)
	{
		build();
	}

	void read(std::ifstream& fin)
	{
		fin>>size_;
		if (size_ == 0) return;
		SizeType tmp;
		fin>>tmp;
		partition_.resize(tmp);
		for (SizeType i = 0; i < tmp; ++i)
			fin>>partition_[i];
		fin>>tmp;
		for (SizeType i = 0; i < tmp; ++i)
			fin>>data_[i];
	}

	void clear()
	{
		partition_.clear();
		data_.clear();
		size_ = 0;
	}

	void resize(SizeType n)
	{
		partition_.clear();
		partition_.resize(2);
		partition_[0] = 0;
		partition_[1] = n;
		size_ = n;
		data_.clear();
		data_.resize(1,0);
	}

	SizeType size() const { return size_; }

	const T& operator[](SizeType i) const
	{
		assert(i < size_);
		SizeType j = findPartition(i);
		return data_[j];
	}

private:

	void build()
	{
		throw PsimagLite::RuntimeError("CompressedVector: build\n");
	}

	SizeType findPartition(SizeType ind) const
	{
		assert(ind < size_);
		for (SizeType i = 1; i < partition_.size(); ++i)
			if (ind < partition_[i]) return i - 1;

		throw PsimagLite::RuntimeError("CompressedVector: findPartition\n");
	}

	SizeType size_;
	PsimagLite::Vector<SizeType>::Type partition_;
	typename PsimagLite::Vector<T>::Type data_;
}; // class CompressedVector

} // namespace Dmrg
#endif // COMPRESSEDVECTOR_H
