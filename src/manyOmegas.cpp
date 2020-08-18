#include <unistd.h>
#include "PsimagLite.h"
#include "Provenance.h"
#include "InputCheck.h"
#include "ManyOmegas.h"
#include "../../dmrgpp/src/Engine/ProgramGlobals.h"

void usage(const PsimagLite::String& name)
{
	std::cerr<<"USAGE is "<<name<<" -f filename [-p precision] [-V]\n";
}

int main(int argc, char** argv)
{
	PsimagLite::PsiApp application("manyOmegas", &argc, &argv, 1);
	typedef Dmrg::ManyOmegas<double> ManyOmegasType;

	int opt = 0;
	bool versionOnly = false;
	PsimagLite::String inputfile;
	PsimagLite::String rootname;
	SizeType precision = 12;
	bool dryrun = false;

	/* PSIDOC DmrgDriver
There is a single input file that is passed as the
argument to \verb!-f!, like so
\begin{lstlisting}
	./dmrg -f input.inp [options]
\end{lstlisting}
The command line arguments
to the main dmrg driver are the following.
	  \begin{itemize}
	  \item[-f] {[}Mandatory, String{]} Dollarized Input to use.
	  \item[-O] {[}Optional, String{]} Root of output file to use.
	  \item[-d] [Optional] Dry run only
	  \item[-p] [Optional, Integer] Digits of precision for printing.
	 \item[-U] [Optional] Make cout output unbuffered
	 \item[-V] [Optional] Print version and exit
	  \end{itemize}
	 */
	while ((opt = getopt(argc, argv,"f:p:O:dV")) != -1) {
		switch (opt) {
		case 'f':
			inputfile = optarg;
			break;
		case 'O':
			rootname = optarg;
			break;
		case 'd':
			dryrun = true;
			break;
		case 'p':
			precision = atoi(optarg);
			std::cout.precision(precision);
			std::cerr.precision(precision);
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

	PsimagLite::String data;
	ManyOmegasType::InputNgType::Writeable::readFile(data, inputfile);
	ManyOmegasType manyOmegas(data, precision, application);

	manyOmegas.run(dryrun, rootname);
}
