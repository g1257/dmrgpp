#include "Concurrency.h"
#include "Engine/BlockDiagonalMatrix.h"
#include "Engine/CmdLineOptions.hh"
#include "Engine/DmrgRunner.h"
#include "Io/IoNg.h"
#include "Provenance.h"

typedef PsimagLite::Vector<PsimagLite::String>::Type VectorStringType;

using namespace Dmrg;

typedef PsimagLite::Concurrency ConcurrencyType;

void printLicense(const PsimagLite::PsiApp& app)
{
	if (!ConcurrencyType::root())
		return;

	std::cout << ProgramGlobals::license;
	Provenance provenance;
	std::cout << provenance;
	std::cout << Provenance::logo(app.name()) << "\n";
	app.checkMicroArch(std::cout, Provenance::compiledMicroArch());
}

void usageOperator()
{
	std::cerr << "USAGE is operator -f filename -e canonical_operator_expression\n";
	std::cerr << "Deprecated options are: -l label [-d dof] [-s site] [-t]\n";
}

int main(int argc, char** argv)
{
	PsimagLite::PsiApp application("DMRG++::dmrg", &argc, &argv, 1);
	PsimagLite::String filename = "";
	int                opt      = 0;
	CmdLineOptions     cmdline_options;
	PsimagLite::String strUsage(application.name());
	strUsage += " -f filename [-k] [-p precision] [-o solverOptions] [-V] [whatToMeasure]";
	bool versionOnly = false;
	/* PSIDOC DmrgDriver
There is a single input file that is passed as the
argument to \verb!-f!, like so
\begin{lstlisting}
	./dmrg -f input.inp whatToMeasure
\end{lstlisting}
where \verb!whatToMeasure! is optional. The command line arguments
to the main dmrg driver are the following.
	  \begin{itemize}
	  \item[-f] {[}Mandatory, String{]} Input to use.
	  \item[-o] {[}Optional, String{]} Extra options for SolverOptions
	  \item[-p] [Optional, Integer] Digits of precision for printing.
	  \item[whatToMeasure] {[}Optional, String{]} What to measure in-situ.
	  This is a comma-separated list of braket specifications.
	  Braket specifications can be bare or dressed, and are explained elsewhere.
	  \item[-l] {[}Optional, String{]} Without this option std::cout is redirected
	  to a file.
	  This option with the string ``?'' prints name of such log file.
	  This option with the string ``-'' writes std::cout to terminal.
	  In other cases, string is the name of the file to redirect std::cout to.
	 \item[-k] [Optional] Keep untar files
	 \item[-U] [Optional] Make cout output unbuffered
	 \item[-S] [Optional, number] Ignore the Threads= line if present in the input,
	  and run with Threads=number
	 \item[-V] [Optional] Print version and exit
	  \end{itemize}
	 */
	while ((opt = getopt(argc, argv, "f:l:p:o:S:kUV")) != -1) {
		switch (opt) {
		case 'f':
			filename = optarg;
			break;
		case 'l':
			cmdline_options.logfile = optarg;
			break;
		case 'p':
			cmdline_options.precision = atoi(optarg);
			std::cout.precision(cmdline_options.precision);
			std::cerr.precision(cmdline_options.precision);
			break;
		case 'o':
			cmdline_options.solver_options += optarg;
			break;
		case 'S':
			cmdline_options.number_of_threads = atoi(optarg);
			break;
		case 'U':
			cmdline_options.unbuffered_output = true;
			break;
		case 'V':
			versionOnly             = true;
			cmdline_options.logfile = "-";
			break;
		default:
			InputCheck::usageMain(strUsage);
			return 1;
		}
	}

	// sanity checks here
	if (filename == "" && !versionOnly) {
		InputCheck::usageMain(strUsage);
		return 1;
	}

	cmdline_options.in_situ_measurements = (optind < argc) ? argv[optind] : "";

	if (cmdline_options.logfile == "?") { // query only
		std::cout << ProgramGlobals::coutName(filename, application.name()) << "\n";
		return 0;
	} else if (cmdline_options.logfile.empty()) {
		cmdline_options.logfile = ProgramGlobals::coutName(filename, application.name());
	}

	// print license
	if (versionOnly) {
		printLicense(application);
		return 0;
	}

	printLicense(application);

	PsimagLite::String data;
	PsimagLite::InputNg<InputCheck>::Writeable::readFile(data, filename);

	DmrgRunner<double> dmrg_runner(application, data, cmdline_options);

	dmrg_runner.doOneRun();
}
