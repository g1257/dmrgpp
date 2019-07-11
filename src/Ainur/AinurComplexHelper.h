#ifndef AINURCOMPLEXHELPER_H
#define AINURCOMPLEXHELPER_H
#include "../Vector.h"

namespace PsimagLite {

template<typename SomeRealType>
class AinurComplexHelper {

public:

	typedef std::complex<SomeRealType> value_type;

	AinurComplexHelper() : t1_(0), t2_(0), counter_(0) {}

	explicit AinurComplexHelper(const SomeRealType& r)
	    : t1_(r), t2_(0), counter_(0) {}

	value_type toComplex() const
	{
		return (counter_ < 2) ? t1_ : value_type(t1_, t2_);
	}

	int end() const { return 0; }

	void insert(int, const SomeRealType& val)
	{
		if (counter_ == 0)
			t1_ = val;
		else if (counter_ == 1)
			t2_ = val;
		if (counter_ > 1)
			err("Wrong complex\n");
		++counter_;
	}

	void insert(int, const value_type& val)
	{
		err("Wrong complex ...\n");
	}

private:

	SomeRealType t1_;
	SomeRealType t2_;
	SizeType counter_;
};

template<typename ComplexOrRealType, bool>
struct MyProxyFor {
	typedef ComplexOrRealType Type;

	static void copy(ComplexOrRealType& dest, const ComplexOrRealType& src)
	{
		dest = src;
	}

	static void copy(std::vector<ComplexOrRealType>& dest,
	                 const std::vector<ComplexOrRealType>& src)
	{
		dest = src;
	}
};

template<typename ComplexOrRealType>
struct MyProxyFor<ComplexOrRealType, true> {

	typedef AinurComplexHelper<typename Real<ComplexOrRealType>::Type> Type;

	static void copy(std::vector<ComplexOrRealType>& dest,
	                 const std::vector<Type>& src)
	{
		dest.clear();
		const SizeType n = src.size();
		if (n == 0) return;
		dest.resize(n);
		for (SizeType i = 0; i < n; ++i)
			dest[i] = src[i].toComplex();
	}

	static void copy(std::vector<ComplexOrRealType>& dest,
	                 const Type& src)
	{
		dest.push_back(src.toComplex());
	}
};

}
#endif // AINURCOMPLEXHELPER_H
