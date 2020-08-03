#include "DmftSolver.h"
#include "Dispersion.h"
#include <unistd.h>
#include "PsimagLite.h"
#include "Provenance.h"
#include "InputCheck.h"
#include "../../dmrgpp/src/Engine/ProgramGlobals.h"

std::streambuf *GlobalCoutBuffer = 0;
std::ofstream GlobalCoutStream;

void restoreCoutBuffer()
{
	if (GlobalCoutBuffer == 0) return;
	GlobalCoutStream.close();
	std::cout.rdbuf(GlobalCoutBuffer);
}

void usage(const PsimagLite::String& name)
{
	std::cerr<<"USAGE is "<<name<<" -f filename [-p precision] [-V]\n";
}

int main(int argc, char** argv)
{
	PsimagLite::PsiApp application("dmft", &argc, &argv, 1);
	typedef PsimagLite::InputNg<Dmft::InputCheck> InputNgType;
	typedef Dmft::DmftSolver<std::complex<double>,  InputNgType> DmftSolverType;
	typedef DmftSolverType::ParamsDmftSolverType ParamsDmftSolverType;
	int opt = 0;
	bool versionOnly = false;
	PsimagLite::String inputfile;
	PsimagLite::String logfile;
	SizeType precision = 12;
	bool unbuffered = false;
	/* PSIDOC DmrgDriver
There is a single input file that is passed as the
argument to \verb!-f!, like so
\begin{lstlisting}
	./dmrg -f input.inp [options]
\end{lstlisting}
The command line arguments
to the main dmrg driver are the following.
	  \begin{itemize}
	  \item[-f] {[}Mandatory, String{]} Input to use.
	  \item[-p] [Optional, Integer] Digits of precision for printing.
	  \item[-l] {[}Optional, String{]} Without this option std::cout is redirected
	  to a file.
	  This option with the string ``?'' prints name of such log file.
	  This option with the string ``-'' writes std::cout to terminal.
	  In other cases, string is the name of the file to redirect std::cout to.
	 \item[-U] [Optional] Make cout output unbuffered
	 \item[-V] [Optional] Print version and exit
	  \end{itemize}
	 */
	while ((opt = getopt(argc, argv,"f:p:l:U:V")) != -1) {
		switch (opt) {
		case 'f':
			inputfile = optarg;
			break;
		case 'p':
			precision = atoi(optarg);
			std::cout.precision(precision);
			std::cerr.precision(precision);
			break;
		case 'l':
			logfile = optarg;
			break;
		case 'U':
			unbuffered = true;
			break;
		case 'V':
			versionOnly = true;
			logfile = "-";
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

	if (logfile != "-") {
		bool queryOnly = (logfile == "?");
		if (logfile == "" || logfile == "?") {
			logfile = Dmrg::ProgramGlobals::coutName(inputfile);
			if (queryOnly) {
				std::cout<<logfile<<"\n";
				return 0;
			}
		}
	}

	if (versionOnly) return 0;

	bool echoInput = false;
	if (logfile != "-") {
		GlobalCoutStream.open(logfile.c_str(), std::ofstream::out);
		if (!GlobalCoutStream || GlobalCoutStream.bad() || !GlobalCoutStream.good()) {
			PsimagLite::String str(application.name());
			str += ": Could not redirect std::cout to " + logfile + "\n";
			err(str);
		}

		echoInput = true;

		std::cerr<<Provenance::logo(application.name());
		std::cerr<<"Standard output sent to ";
		std::cerr<<logfile<<"\n";
		std::cerr.flush();
		GlobalCoutBuffer = std::cout.rdbuf(); //save old buf
		std::cout.rdbuf(GlobalCoutStream.rdbuf()); //redirect std::cout to file
		if (unbuffered) {
			std::cout.setf(std::ios::unitbuf);
			GlobalCoutStream.setf(std::ios::unitbuf);
		}

		atexit(restoreCoutBuffer);
	}

	application.printCmdLine(std::cout);
	if (echoInput) application.echoBase64(std::cout, inputfile);

	Dmft::InputCheck inputCheck;
	InputNgType::Writeable ioWriteable(inputfile, inputCheck);
	InputNgType::Readable io(ioWriteable);

	ParamsDmftSolverType params(io);
	if (precision > 0) params.precision = precision;
	params.echoInput = echoInput;

	DmftSolverType dmftSolver(params, application);

	dmftSolver.selfConsistencyLoop();
}
