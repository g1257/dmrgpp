#include <iostream>
#include <fstream>
#include "Vector.h"

void procFile(PsimagLite::String filename)
{
	std::ifstream fin(filename.c_str());
	// replace ""
	// replace ''
	// replace top level {}

}

int main(int argc, char** argv)
{
	if (argc == 1) return 1;
	procFile(argv[1]);
}
