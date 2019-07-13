#include "PredicateAwesome.h"
#include <cstdlib>

class BogusClass {};

int main(int argc, char** argv)
{
	if (argc < 3) return 1;
	PsimagLite::String predicate(argv[1]);

	PsimagLite::PredicateAwesome<BogusClass> pAwesome(predicate);
	std::cout<<pAwesome.isTrue("l", atoi(argv[2]));
}
