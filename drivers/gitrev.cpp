#include "GitRevision.h"

int main(int,char **)
{
	PsimagLite::GitRevision gitrev("./","PsimagLite");
	std::cout<<gitrev;

}
