#ifndef PSI_OPTIONAL_H
#define PSI_OPTIONAL_H

#include <cassert>

namespace PsimagLite {
template<typename T>
class Optional {

public:

	Optional(T& t) : t_(&t)
	{}

	Optional(int) : t_(0)
	{}

	bool nonNull() const
	{
		return (t_);
	}

	T& operator()()
	{
		assert(t_);
		return *t_;
	}

private:

	T* t_;
};
}

#endif // PSI_OPTIONAL_H
