#ifndef TEST_MEM_RESOLV_1_H
#define TEST_MEM_RESOLV_1_H
#include <iostream>
#include <cstdlib>
#include "MemResolv.h"
#include "IsClass.h"

typedef PsimagLite::MemResolv MemResolv;

template<typename T>
class TestMemResolv1 {

	static const bool IS_CLASS = PsimagLite::IsClass<T>::value;
	typedef std::vector<T> VectorType;
	typedef PsimagLite::ResolveFinalOrNot<VectorType,IS_CLASS> ResolveFinalOrNotType;
	typedef TestMemResolv1<T> ThisType;

public:

	TestMemResolv1(int size)
	: size_(size),data_(size)
	{}

	void setTo(T x)
	{
		for (int i = 0; i < data_.size(); ++i)
			data_[i] = x;
	}

	const T& get(int i) const
	{
		return data_[i];
	}

	SizeType memResolv(MemResolv& vmptr,
	                   SizeType x = 0,
	                   PsimagLite::String msg = "") const
	{
		PsimagLite::String str = msg;
		msg += "TestMemResolv1 ";
		const char* start = (const char*)&size_;
		const char* end = (const char*)&data_;
		SizeType total = vmptr.memResolv(&size_,end - start,str + "size");
		total += vmptr.memResolv(&data_,size_,str + "size");
		return total;
	}

	template<typename T2>
	friend std::ostream& operator<<(std::ostream& os, const TestMemResolv1<T2>& mt);

private:

	TestMemResolv1(const TestMemResolv1& other);

	const TestMemResolv1& operator=(const TestMemResolv1& other);

	int size_;
	VectorType data_;
};

template<typename T>
std::ostream& operator<<(std::ostream& os, const TestMemResolv1<T>& mt)
{
	os<<"TestMemResolv1 size= "<<mt.size_<<"\n";
	for (SizeType i = 0; i < mt.size_; ++i)
		os<<mt.data_[i]<<" ";
	os<<"\n";
	return os;
}

#endif

