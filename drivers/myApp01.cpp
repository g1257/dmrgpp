#include "PsimagLite.h"

using namespace PsimagLite;

// g++ myApp01.cpp -o myApp01 -I../src -L../lib -lpsimaglite
int main(int argc, char** argv)
{
	SizeType nthreads = 1;
	PsiApp application("DMRG++",&argc,&argv,nthreads);
}
