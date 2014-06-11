#ifndef PSI_IS_CLASS_H
#define PSI_IS_CLASS_H
#include<iostream>

namespace PsimagLite {

template<typename T>
class IsClass {

	typedef char One;
	typedef struct { char a[2]; } Two;

	template<typename C>
	static One test(int C::*);

	template<typename C>
	static Two test(...);

public:

	enum { value = sizeof(IsClass<T>::template test<T>(0)) == sizeof(One) };

};

}

#endif

