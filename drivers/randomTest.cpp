#include "Random48.h"
#include <cstdlib>
#include <iostream>
using namespace PsimagLite;

typedef double RealType;
typedef Random48<RealType> RandomType;

int main(int, char* argv[])
{
	RandomType rng(100);
	rng.seed(std::atoi(argv[1]));
	RealType x = rng.random();
	std::cout << "x=" << x << "\n";
}
