#include "CincuentaInputCheck.h"
#include "Dispersion.h"
#include "DmftSolver.h"
#include "ImpuritySolverNeqTdmrg.h"
#include "NeqDmftSolver.h"
#include "ProgramGlobals.h"
#include "Provenance.h"
#include <PsimagLite/InputPath.hpp>
#include <PsimagLite/PsimagLite.h>
#include <unistd.h>

std::streambuf* GlobalCoutBuffer = 0;
std::ofstream   GlobalCoutStream;

void restoreCoutBuffer()
{
	if (GlobalCoutBuffer == 0)
		return;
	GlobalCoutStream.close();
	std::cout.rdbuf(GlobalCoutBuffer);
}

void usage(const std::string& name)
{
	std::cerr << "USAGE is " << name << " -f filename [-p precision] [-V]\n";
}

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
  \item[-I] {Optional, String} Add an input path to the search; it can be used multiple times
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
int main(int argc, char** argv)
{
	PsimagLite::PsiApp application("dmft", &argc, &argv, 1);
	using InputNgType                 = PsimagLite::InputNg<Dmft::CincuentaInputCheck>;
	using RealType                    = double;
	using DmftSolverType              = Dmft::DmftSolver<std::complex<RealType>>;
	using ParamsDmftSolverType        = DmftSolverType::ParamsDmftSolverType;
	int                   opt         = 0;
	bool                  versionOnly = false;
	std::string           inputfile;
	std::string           logfile;
	SizeType              precision  = 12;
	bool                  unbuffered = false;
	PsimagLite::InputPath input_path;

	while ((opt = getopt(argc, argv, "f:p:l:U:I:V")) != -1) {
		switch (opt) {
		case 'f':
			inputfile = optarg;
			break;
		case 'I':
			input_path.push(optarg);
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
			logfile     = "-";
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

	using ConcurrencyType = PsimagLite::Concurrency;

	// print license
	if (ConcurrencyType::root()) {
		Provenance provenance;
		std::cout << provenance;
		std::cout << Provenance::logo(application.name()) << "\n";
		application.checkMicroArch(std::cout, Provenance::compiledMicroArch());
	}

	if (logfile != "-") {
		bool queryOnly = (logfile == "?");
		if (logfile == "" || logfile == "?") {
			logfile = Dmrg::ProgramGlobals::coutName(inputfile, "cincuenta");
			if (queryOnly) {
				std::cout << logfile << "\n";
				return 0;
			}
		}
	}

	if (versionOnly)
		return 0;

	bool echoInput = false;
	if (logfile != "-") {
		GlobalCoutStream.open(logfile.c_str(), std::ofstream::out);
		if (!GlobalCoutStream || GlobalCoutStream.bad() || !GlobalCoutStream.good()) {
			std::string str(application.name());
			str += ": Could not redirect std::cout to " + logfile + "\n";
			err(str);
		}

		echoInput = true;

		std::cerr << Provenance::logo(application.name());
		std::cerr << "Standard output sent to ";
		std::cerr << logfile << "\n";
		std::cerr.flush();
		GlobalCoutBuffer = std::cout.rdbuf(); // save old buf
		std::cout.rdbuf(GlobalCoutStream.rdbuf()); // redirect std::cout to file
		if (unbuffered) {
			std::cout.setf(std::ios::unitbuf);
			GlobalCoutStream.setf(std::ios::unitbuf);
		}

		atexit(restoreCoutBuffer);
	}

	application.printCmdLine(std::cout);
	if (echoInput)
		application.echoBase64(std::cout, inputfile);

	Dmft::CincuentaInputCheck inputCheck;
	InputNgType::Writeable    ioWriteable(input_path.findFirst(inputfile), inputCheck);
	InputNgType::Readable     io(ioWriteable);

	ParamsDmftSolverType params(io);
	// BEGIN adjust params
	if (precision > 0) {
		params.precision = precision;
	}

	params.echoInput = echoInput;

	// END adjust params
	DmftSolverType::FitType::InitResults initResults(io);

	DmftSolverType dmftSolver(params, initResults, application, io);

	dmftSolver.selfConsistencyLoop();

	dmftSolver.print(std::cout);

	// Non-equilibrium DMFT mode: triggered when TmaxNeq is present in input.
	// The equilibrium DMFT run above supplies the bath parameters (fixed bath
	// approximation for the interaction quench U_i -> U_f).
	{
		using ParamsNeqType = Dmft::ParamsNeqDmftSolver<std::complex<RealType>>;

		RealType tMax  = 0;
		bool     isNeq = false;
		try {
			io.readline(tMax, "TmaxNeq=");
			isNeq = (tMax > 0);
		} catch (std::exception&) { }

		if (isNeq) {
			std::cout << "\n=== Non-equilibrium DMFT (interaction quench) ===\n";
			ParamsNeqType neqParams(io);

			// Check for tDMRG solver selection
			std::string neqSolverType;
			try {
				io.readline(neqSolverType, "NeqSolver=");
			} catch (std::exception&) { }

			if (neqSolverType == "tdmrg") {
				std::cout << "  using ImpuritySolverNeqTdmrg (tDMRG)\n";
				using TdmrgImpType
				    = Dmft::ImpuritySolverNeqTdmrg<std::complex<RealType>>;
				TdmrgImpType tdmrgSolver(neqParams, application, io);
				tdmrgSolver.solve(dmftSolver.bathResult());
				const std::string& p = neqParams.neqOutputPrefix;
				tdmrgSolver.gimp().dump(p.empty() ? "green" : p + "-green");
			} else {
				using ExactNeqSolverType
				    = Dmft::NeqDmftSolver<std::complex<RealType>>;
				ExactNeqSolverType neqSolver(neqParams, io);
				neqSolver.solve(dmftSolver.bathResult());
				neqSolver.dumpGreenFunctions();
			}
		}
	}
}
