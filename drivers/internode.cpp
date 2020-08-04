#include "InterNode.h"

int main(int argc, char* argv[])
{
	if (argc < 2) {
		std::cerr<<"USAGE "<<argv[0]<<" number\n";
		return 1;
	}

	const SizeType n = atoi(argv[1]);

	// original
	for (SizeType i = 0; i < n; ++i) {
		std::cout<<i;
	}

	std::cout<<"\n--------------------------\n";

	//lambda
	PsimagLite::InterNode<> internode;

	internode.parallelFor(0, n, [](SizeType i, SizeType){std::cout << i;});
	std::cout<<"\n--------------------------\n";

}
