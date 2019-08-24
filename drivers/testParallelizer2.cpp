#include "Parallelizer2.h"
#include "Vector.h"
#include "LoadBalancerWeights.h"

int main(int argc, char** argv)
{
	if (argc < 3) return 1;

	const SizeType n = atoi(argv[1]);
	const SizeType threads = atoi(argv[2]);

	PsimagLite::Vector<double>::Type v(n);

	PsimagLite::Parallelizer2<> parallelizer(threads);

	std::cout<<"Testing Parallelizer2 with "<<parallelizer.name();
	std::cout<<" and "<<parallelizer.numberOfThreads()<<" threads.\n";
	parallelizer.parallelFor(0, n, [&v](SizeType i, SizeType)
	{
		v[i] = i + 42; 
	});

	/*
	for (SizeType i = 0; i < n; ++i) {
		v[i] = i + 42;
	}
	*/

	std::cout<<v[0]<<" "<<v[n - 1]<<"\n";
}
