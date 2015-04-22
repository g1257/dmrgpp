#include "GitRevision.h"

int main(int,char **)
{
	PsimagLite::GitRevision gitrev("./","dmrgpp");
	std::cout<<gitrev;
	PsimagLite::GitRevision gitrev2("../../PsimagLite/","psimagLite");
	std::cout<<gitrev2;

}
