#include "Ainur.h"

void partiallyReadSomething(const PsimagLite::Ainur& ainur)
{
	SizeType n = 0;
	ainur.readValue(n, "TotalNumberOfSites");
	std::cout<<"TotalNumberOfSites="<<n<<"\n";
}

int main(int argc, char** argv)
{
	if (argc == 1) return 1;
	std::ifstream fin(argv[1]);
	PsimagLite::String str;

	fin.seekg(0, std::ios::end);
	str.reserve(fin.tellg());
	fin.seekg(0, std::ios::beg);

	str.assign((std::istreambuf_iterator<char>(fin)),
	           std::istreambuf_iterator<char>());
	fin.close();

	PsimagLite::Ainur ainur(str);
	partiallyReadSomething(ainur);
	ainur.printAll(std::cout);
}
