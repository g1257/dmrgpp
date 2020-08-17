#include <unistd.h>
#include "PsimagLite.h"
#include "Provenance.h"
#include "InputCheck.h"
#include "ProcOmegas.h"
#include "../../dmrgpp/src/Engine/ProgramGlobals.h"

void usage(const PsimagLite::String& name)
{
	std::cerr<<"USAGE is "<<name<<" -f filename [other options]\n";
	std::cerr<<"-f filename; Mandatory template filename\n";
	std::cerr<<"-p precision; Precision in number of decimals\n";
	std::cerr<<"-X; Skip Fourier transform and only gather data\n";
	std::cerr<<"-I {[}Optional, String{]} Root of input file to use.\n";
	std::cerr<<"-O {[}Optional, String{]} Root of output file to use.\n";
	std::cerr<<"-V; Print version and exit\n";
	std::cerr<<"\nLimitations: Ainur only; No Cheby yet\n";
}

int main(int argc, char** argv)
{
	PsimagLite::PsiApp application("procOmegas", &argc, &argv, 1);
	typedef Dmrg::ProcOmegas<double> ProcOmegasType;
	typedef ProcOmegasType::InputNgType InputNgType;

	int opt = 0;
	bool versionOnly = false;
	PsimagLite::String inputfile;
	PsimagLite::String rootIname = "input";
	PsimagLite::String rootOname = "out";
	SizeType precision = 12;
	bool skipFourier = false;

	/* PSIDOC DmrgDriver
There is a single input file that is passed as the
argument to \verb!-f!, like so
\begin{lstlisting}
	./procOmegas -f input.inp [options]
\end{lstlisting}
The command line arguments
to the main dmrg driver are the following.
	  \begin{itemize}
	  \item[-f] {[}Mandatory, String{]} Dollarized Input to use.
	  \item[-I] {[}Optional, String{]} Root of input file to use.
	  \item[-O] {[}Optional, String{]} Root of output file to use.
	  \item[-X] [Optional] Skip Fourier transform and only gather data
	  \item[-p] [Optional, Integer] Digits of precision for printing.
	 \item[-V] [Optional] Print version and exit
	  \end{itemize}
	 */
	while ((opt = getopt(argc, argv,"f:p:I:O:XV")) != -1) {
		switch (opt) {
		case 'f':
			inputfile = optarg;
			break;
		case 'I':
			rootIname = optarg;
			break;
		case 'O':
			rootOname = optarg;
			break;
		case 'X':
			skipFourier = true;
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

	Dmrg::InputCheck inputCheck;
	InputNgType::Writeable ioW(inputfile, inputCheck);
	InputNgType::Readable io(ioW);
	ProcOmegasType procOmegas(io,
	                          precision,
	                          skipFourier,
	                          rootIname,
	                          rootOname);

	procOmegas.run();

	procOmegas.printPgfplots(rootOname + ".pgfplots");
}
