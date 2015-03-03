#include "SumDecomposition.h"

int main(int argc, char* argv[])
{
	if (argc < 3) return 1;
	int total = atoi(argv[1]);
	int sum = atoi(argv[2]);
	int selection = (argc > 3) ? atoi(argv[3]) : 0;
	PsimagLite::SumDecomposition::SelEnum sel = PsimagLite::SumDecomposition::SEL_ALL;
	if (selection < 0) sel = PsimagLite::SumDecomposition::SEL_SIZE;
	if (selection >=0 && argc > 3) sel = PsimagLite::SumDecomposition::SEL_INDEX;
	PsimagLite::SumDecomposition sd(total, sum, sel, selection);
	if (sel != PsimagLite::SumDecomposition::SEL_SIZE)
		std::cout<<sd;
	else
		std::cout<<sd.size()<<"\n";
}

