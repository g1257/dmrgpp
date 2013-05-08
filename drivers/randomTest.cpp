#include <iostream>
#include <cstdlib>
#include "Random48.h"
using namespace PsimagLite;

typedef double RealType;
typedef Random48<RealType> RandomType;

int main(int argc,char* argv[])
{
	RandomType rng(100);
	rng.seed(std::atoi(argv[1]));
	RealType x = rng.random();
	std::cout<<"x="<<x<<"\n";
}

