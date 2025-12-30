#include "Concurrency.h"
#include "Random48.h"
#include <cstdlib>
#include <iostream>

using namespace PsimagLite;

typedef double             RealType;
typedef Random48<RealType> RandomType;

int main(int argc, char* argv[])
{
	constexpr unsigned int nthreads = 1;
	PsimagLite::Concurrency(&argc, &argv, nthreads);

	if (argc != 2) {
		std::cerr << "USAGE: " << argv[0] << " seed\n";
		return 1;
	}

	RandomType rng(100);
	rng.seed(std::atoi(argv[1]));
	RealType x = rng.random();
	std::cout << "x=" << x << "\n";
}
