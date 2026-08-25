#include "CincuentaInputCheck.h"
#include "Dispersion.h"
#include "DmftSolver.h"
#include "ImpuritySolverNeqGBEK.h"
#include "ImpuritySolverNeqTdmrg.h"
#include "NeqDmftSolver.h"
#include "ProgramGlobals.h"
#include "Provenance.h"
#include <PsimagLite/InputPath.hpp>
#include <PsimagLite/PsimagLite.h>
#include <unistd.h>

class NullStreamBuffer : public std::streambuf {
protected:

	int overflow(int c) override { return c; }
};

std::streambuf*  GlobalCoutBuffer = 0;
std::ofstream    GlobalCoutStream;
NullStreamBuffer GlobalNullStream;

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

// Reject presence of an equilibrium-bath-fit-only key under NeqAtomicLimit=1,
// where the equilibrium DMFT stage never runs at all, so the key would
// otherwise be silently parsed and ignored. T must match the key's declared
// Ainur type (see CincuentaInputCheck::import()).
template <typename T, typename ReadableType>
void rejectUnderAtomicLimit(ReadableType& io, const std::string& key)
{
	T    tmp {};
	bool present = true;
	try {
		io.readline(tmp, key);
	} catch (std::exception&) {
		present = false;
	}
	if (present)
		err(key
		    + " has no effect when NeqAtomicLimit=1 (the equilibrium bath "
		      "fit is skipped entirely); remove it from the input\n");
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
	const bool isRoot     = ConcurrencyType::root();

	// print license
	if (isRoot) {
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
				if (isRoot)
					std::cout << logfile << "\n";
				return 0;
			}
		}
	}

	if (versionOnly)
		return 0;

	bool echoInput = false;
	if (!isRoot) {
		// Only rank zero owns user-facing output and shared log files.
		GlobalCoutBuffer = std::cout.rdbuf();
		std::cout.rdbuf(&GlobalNullStream);
		atexit(restoreCoutBuffer);
	} else if (logfile != "-") {
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

	if (isRoot) {
		application.printCmdLine(std::cout);
		if (echoInput)
			application.echoBase64(std::cout, inputfile);
	}

	Dmft::CincuentaInputCheck inputCheck;
	InputNgType::Writeable    ioWriteable(input_path.findFirst(inputfile), inputCheck);
	InputNgType::Readable     io(ioWriteable);

	// NeqAtomicLimit=1 means there is no equilibrium DMFT stage at all: the
	// neq run starts from an empty bath instead of a fitted one. Peek it
	// before touching anything equilibrium-specific, so that case can skip
	// *constructing* (not just running) the equilibrium machinery, and
	// reject equilibrium-only keywords that would otherwise silently do
	// nothing (see cincuenta/doc/neq_dmft_ed_input.md).
	bool neqAtomicLimit = false;
	{
		int tmp = 0;
		try {
			io.readline(tmp, "NeqAtomicLimit=");
			neqAtomicLimit = (tmp > 0);
		} catch (std::exception&) { }
	}

	PsimagLite::Vector<RealType>::Type equilibriumBathResult;

	if (neqAtomicLimit) {
		// These configure the equilibrium bath fit, which never runs in this
		// mode -- reject them outright rather than silently ignoring them.
		// (LatticeGf= is still required, not rejected: it's checked inside
		// ParamsMatsubaraGrid, constructed below as part of ParamsNeqDmftSolver,
		// which only accepts the zero-bandwidth value under NeqAtomicLimit=1.)
		rejectUnderAtomicLimit<SizeType>(io, "NumberOfBathPoints=");
		rejectUnderAtomicLimit<std::string>(io, "ImpuritySolver=");
		rejectUnderAtomicLimit<std::string>(io, "FitOptions=");

		std::cout << "\nNeqAtomicLimit=1: no equilibrium DMFT stage -- the neq "
		             "run starts from an empty bath.\n";
	} else {
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

		if (isRoot)
			dmftSolver.print(std::cout);

		equilibriumBathResult = dmftSolver.bathResult();
	}

	// Non-equilibrium DMFT mode: triggered when TmaxNeq is present in input.
	// The equilibrium DMFT stage above supplies the bath parameters (fixed
	// bath approximation for the interaction quench U_i -> U_f), unless
	// NeqAtomicLimit=1, in which case equilibriumBathResult is empty.
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

			// NeqSolver= selects the neq propagation method: "exactdiag"
			// (default, also selected by omitting NeqSolver= entirely) or
			// "tdmrg".
			// "exactdiag" always uses the GBEK-capable exact-diagonalization solver:
			// NeqBathRank=0 (default) reduces exactly to a fixed equilibrium
			// bath with no time-dependent second bath; NeqBathRank>0 adds the
			// low-rank Cholesky second bath (Gramsch/Balzer/Eckstein/Kollar,
			// PRB 88, 235106 (2013)) so the bath can respond to the quench.
			std::string neqSolverType;
			try {
				io.readline(neqSolverType, "NeqSolver=");
			} catch (std::exception&) { }

			if (neqSolverType == "tdmrg") {
				if (neqParams.neqAtomicLimit)
					err("NeqSolver=\"tdmrg\" does not support NeqAtomicLimit=1 "
					    "(untested combination)\n");
				std::cout << "  using ImpuritySolverNeqTdmrg (tDMRG)\n";
				using TdmrgImpType
				    = Dmft::ImpuritySolverNeqTdmrg<std::complex<RealType>>;
				TdmrgImpType tdmrgSolver(neqParams, application, io);
				tdmrgSolver.solve(equilibriumBathResult);
				const std::string& p = neqParams.neqOutputPrefix;
				tdmrgSolver.gimp().dump(p.empty() ? "green" : p + "-green");
			} else if (neqSolverType == "" || neqSolverType == "exactdiag") {
				std::cout << "  using ImpuritySolverNeqGBEK (exact "
				             "diagonalization, NeqBathRank="
				          << neqParams.neqBathRank << ")\n";
				using EdNeqSolverType
				    = Dmft::NeqDmftSolver<std::complex<RealType>,
				                          Dmft::ImpuritySolverNeqGBEK>;
				EdNeqSolverType neqSolver(neqParams, io);
				neqSolver.solve(equilibriumBathResult);
				neqSolver.dumpGreenFunctions();
			} else {
				err("Unknown NeqSolver=\"" + neqSolverType
				    + "\"; expected \"exactdiag\" (default) or \"tdmrg\"\n");
			}
		}
	}
}
