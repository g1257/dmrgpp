#include "Ainur.h"

PsimagLite::String dmrgImport()
{
	PsimagLite::String str("");
	str += "require TotalNumberOfSites;";
	return str;
}

int main(int argc, char** argv)
{
	if (argc == 1) return 1;
	PsimagLite::String dmrgppImport = dmrgImport();
	PsimagLite::Ainur ainur(argv[1], dmrgppImport);
}
