#include <iostream>
#include "TarPack.h"
#include "InputNg.h"
#include "InputCheck.h"
#include "ParametersDmrgSolver.h"
#include "ToolBox.h"

#ifndef USE_FLOAT
typedef double RealType;
#else
typedef float RealType;
#endif
typedef PsimagLite::InputNg<Dmrg::InputCheck> InputNgType;
typedef Dmrg::ParametersDmrgSolver<RealType,InputNgType::Readable>
ParametersDmrgSolverType;

void usage(const char* name)
{
	std::cerr<<"USAGE is "<<name<<" -f filename -a action [-p precision]\n";
}

int main(int argc,char *argv[])
{
	using namespace Dmrg;

	PsimagLite::String filename;
	PsimagLite::String action;
	int opt = 0;
	int precision = 6;
	while ((opt = getopt(argc, argv,"f:p:a:")) != -1) {
		switch (opt) {
		case 'f':
			filename = optarg;
			break;
		case 'p':
			precision = atoi(optarg);
			std::cout.precision(precision);
			std::cerr.precision(precision);
			break;
		case 'a':
			action = optarg;
			break;
		default:
			usage(argv[0]);
			return 1;
		}
	}

	PsimagLite::String list = (optind < argc) ? argv[optind] : "";

	//sanity checks here
	if (filename=="" || action == "") {
		usage(argv[0]);
		return 1;
	}

	typedef PsimagLite::Concurrency ConcurrencyType;
	ConcurrencyType concurrency(&argc,&argv,1);

	// print license
	if (ConcurrencyType::root()) {
		std::cerr<<ProgramGlobals::license;
		Provenance provenance;
		std::cout<<provenance;
	}

	InputCheck inputCheck;
	InputNgType::Writeable ioWriteable(filename,inputCheck);
	InputNgType::Readable io(ioWriteable);

	//! Read the parameters for this run
	bool earlyExit = true;
	ParametersDmrgSolverType dmrgSolverParams(io,earlyExit);

	ConcurrencyType::npthreads = dmrgSolverParams.nthreads;

	if (Dmrg::ToolBox::actionCanonical(action) == Dmrg::ToolBox::ACTION_ENERGIES) {
		Dmrg::ToolBox::printEnergies(dmrgSolverParams.filename);
	} else {
		std::cerr<<argv[0]<<": Unknown action "<<action<<"\n";
		std::cerr<<"\tSupported actions are "<<Dmrg::ToolBox::actions()<<"\n";
		return 1;
	}
} // main

