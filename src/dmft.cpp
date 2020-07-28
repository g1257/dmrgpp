#include "DmftSolver.h"
#include "Dispersion.h"
#include <unistd.h>
#include "PsimagLite.h"
#include "Provenance.h"
#include "InputCheck.h"

void usage(const PsimagLite::String& name)
{
	std::cerr<<"USAGE is "<<name<<" -f filename [-p precision] [-V]\n";
}

int main(int argc, char** argv)
{
	PsimagLite::PsiApp application("observe",&argc,&argv,1);
	typedef Dmft::DmftSolver<std::complex<double> > DmftSolverType;
	typedef DmftSolverType::RealType RealType;
	typedef DmftSolverType::DispersionType DispersionType;
	typedef PsimagLite::InputNg<Dmft::InputCheck> InputNgType;

	int opt = 0;
	bool versionOnly = false;
	PsimagLite::String inputfile;

	while ((opt = getopt(argc, argv,"f:V")) != -1) {
		switch (opt) {
		case 'f':
			inputfile = optarg;
			break;
		case 'V':
			versionOnly = true;
			break;
		default:
			usage(application.name());
			return 1;
		}
	}

	if (inputfile == "") {
		if (!versionOnly) {
			usage(application.name());
			return 1;
		}
	}

	typedef PsimagLite::Concurrency ConcurrencyType;

	// print license
	if (ConcurrencyType::root()) {
		Provenance provenance;
		std::cout<<provenance;
		std::cout<<Provenance::logo(application.name())<<"\n";
		application.checkMicroArch(std::cout, Provenance::compiledMicroArch());
	}

	if (versionOnly) return 0;

	application.printCmdLine(std::cout);
	Dmft::InputCheck inputCheck;
	InputNgType::Writeable ioWriteable(inputfile, inputCheck);
	InputNgType::Readable io(ioWriteable);

	DispersionType dispersion(numberOfKpoints);

	DmftSolverType dmftSolver(ficticiousBeta, nMatsubara, dispersion, mu);

	dmftSolver.selfConsistencyLoop();
}
