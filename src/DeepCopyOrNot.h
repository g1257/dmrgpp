#ifndef DEEPCOPYORNOT_H
#define DEEPCOPYORNOT_H
#include "Vector.h"

namespace PsimagLite {

template<typename UnderlyingType> // UnderlyingType must have copy ctor
class DeepCopyOrNot {

public:

	typedef typename Vector<UnderlyingType*>::Type VectorUnderlyingType;

	DeepCopyOrNot(bool isDeep) : isDeep_(isDeep) {}

	~DeepCopyOrNot()
	{
		SizeType n = garbage_.size();
		for (SizeType i = 0; i < n; ++i) {
			delete garbage_[i];
			garbage_[i] = 0;
		}
	}

	const UnderlyingType& operator()(const UnderlyingType& other)
	{
		if (!isDeep_) return other;

		UnderlyingType* tmp = new UnderlyingType(other);
		garbage_.push_back(tmp);
		return *tmp;
	}

private:

	DeepCopyOrNot(const DeepCopyOrNot&);

	DeepCopyOrNot& operator=(const DeepCopyOrNot&);

	bool isDeep_;
	VectorUnderlyingType garbage_;
};
}
#endif // DEEPCOPYORNOT_H
