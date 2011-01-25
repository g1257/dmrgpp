#include <iostream>
#include <cstdlib>
#include "Random48.h"
using namespace PsimagLite;

typedef double RealType;
typedef Random48<RealType> RandomType;

int main(int argc,char* argv[])
{
	RandomType::seed(std::atoi(argv[1]));
	RealType x = RandomType::random();
	std::cout<<"x="<<x<<"\n";
}

