#ifndef ARRAY_H
#define ARRAY_H
#include "AllocatorCpu.h"

namespace Dmrg {

template<typename T>
class Array {

public:

	Array() : size_(0), data_(0) {}

	Array(SizeType n) : size_(0), data_(0)
	{
		allocate(n);
	}

	Array(const Array& other) : size_(0), data_(0)
	{
		clear();
		allocate(other.size_);
		memcpy(data_, other.data_, size_*sizeof(T));
	}

	~Array()
	{
		delete [] data_;
		data_ = 0;
	}


	Array& operator=(const Array& other)
	{
		clear();
		allocate(other.size_);
		memcpy(data_, other.data_, size_*sizeof(T));
		return *this;
	}

	void clear()
	{
		size_ = 0;
		delete [] data_;
		data_ = 0;
	}

	void resize(SizeType n)
	{
		if (size_ == n) return;
		clear();
		allocate(n);
	}

	SizeType size() const { return size_; }

	const T& operator[](SizeType ind) const
	{
		assert(ind < size_);
		return data_[ind];
	}

	T& operator[](SizeType ind)
	{
		assert(ind < size_);
		return data_[ind];
	}

private:

	void allocate(SizeType n)
	{
		assert(data_ == 0 && size_ == 0);
		size_ = n;
		data_ = new T[size_];
	}

	SizeType size_;
	T* data_;
};
}
#endif // ARRAY_H
