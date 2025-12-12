#include "Concurrency.h"
#include "PredicateAwesome.h"
#include <cstdlib>

int main(int argc, char** argv)
{
	constexpr unsigned int nthreads = 1;
	PsimagLite::Concurrency(&argc, &argv, nthreads);

	if (argc != 3) {
		std::cerr << "USAGE: " << argv[0] << " predicate value_of_l\n";
		return 1;
	}

	PsimagLite::String predicate(argv[1]);

	PsimagLite::replaceAll(predicate, "c", "5");
	PsimagLite::PredicateAwesome<> pAwesome(predicate);
	std::cout << pAwesome.isTrue("l", atoi(argv[2]));
}
