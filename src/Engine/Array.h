#ifndef ARRAY_H
#define ARRAY_H
#include "AllocatorCpu.h"
#include "Vector.h"
#include <cstring>

namespace Dmrg {

template<typename T>
class Array {

public:

	typedef T value_type;

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

	Array(const std::vector<T>& other) : size_(0), data_(0)
	{
		fromStdVector(other);
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

	void toStdVector(std::vector<T>& v) const
	{
		if (size_ == 0) return;
		v.resize(size_);
		memcpy(&(v[0]), data_, sizeof(SizeType)*size_);
	}

	void fromStdVector(const std::vector<T>& v)
	{
		clear();
		allocate(v.size());
		memcpy(data_, &(v[0]), size_*sizeof(T));
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
