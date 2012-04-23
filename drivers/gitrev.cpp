#include "GitRevision.h"

int main(int argc,char *argv[])
{
	PsimagLite::GitRevision gitrev("./","PsimagLite");
	std::cout<<gitrev;

}
