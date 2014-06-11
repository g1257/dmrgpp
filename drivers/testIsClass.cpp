#include "IsClass.h"
#include <vector>

class A {};

template<typename T>
int testIsClass()
{
	return PsimagLite::IsClass<T>::True;
}

int main()
{
	std::cout<<testIsClass<A>()<<"\n";
	std::cout<<testIsClass<double>()<<"\n";
	std::cout<<testIsClass<int>()<<"\n";
	std::cout<<testIsClass<std::vector<double> >()<<"\n";
	std::cout<<testIsClass<char>()<<"\n";
}

