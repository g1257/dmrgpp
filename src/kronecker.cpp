#include <iostream>
#include <unistd.h> // for getopt
#include "Vector.h"
#include "Concurrency.h"
#include "ProgramGlobals.h"
#include "Provenance.h"
#include "MiniAppKronecker.h"

void usage(const char* name)
{
	std::cerr<<"USAGE is "<<name<<" -f filename | -V\n";
}

int main(int argc, char** argv)
{
	int opt = 0;
	PsimagLite::String filename("");
	bool versionOnly = false;
	while ((opt = getopt(argc, argv,"f:V")) != -1) {
		switch (opt) {
		case 'f':
			filename = optarg;
			break;
		case 'V':
			versionOnly = true;
			break;
		default:
			usage(argv[0]);
			return 1;
		}
	}

	//sanity checks here
	if (filename == "") {
		if (!versionOnly) {
			usage(argv[0]);
			return 1;
		}
	}

	typedef PsimagLite::Concurrency ConcurrencyType;
	ConcurrencyType concurrency(&argc,&argv,1);

	// print license
	if (ConcurrencyType::root()) {
		std::cerr<<Dmrg::ProgramGlobals::license;
		Provenance provenance;
		std::cout<<provenance;
	}

	if (versionOnly) return 0;
	Dmrg::MiniAppKronecker miniAppKron(filename);
	miniAppKron.check();
}
