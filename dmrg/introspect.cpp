#include "Concurrency.h"
#include "Engine/BlockDiagonalMatrix.h"
#include "Engine/CmdLineOptions.hh"
#include "Engine/DmrgRunner.h"
#include "Engine/OptionsForIntrospect.hh"
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

int main(int argc, char** argv)
{
	PsimagLite::PsiApp   application("DMRG++::instrospect", &argc, &argv, 1);
	PsimagLite::String   filename = "";
	int                  opt      = 0;
	OptionsForIntrospect operator_options;
	CmdLineOptions       cmdline_options;
	cmdline_options.solver_options = ",introspect";
	PsimagLite::String strUsage(application.name());
	strUsage += " -f filename [-p precision] [-s site] [-V] -e expression [-H | -B]";
	bool versionOnly = false;
	/* PSIDOC OperatorDriver
	 The arguments to the \verb!operator! executable are as follows.
	\begin{itemize}
	 \item[-f] [Mandatory, String] Input to use. The Model= line is
	very important in input.inp.

	\item[-e] [Mandatory unless -H or -B, String] OperatorExpression; see manual

	\item[-s] [Optional, Integer] \emph{Deprecated. Use -e.}
	Site for operator.
	Meaningful only for Models where
	the Hilbert space depends on the site (different kinds of atoms).
	Defaults to 0.

	\item[-B] [Optional] Prints the basis and all operators for the model

	\item[-H] [Optional] Prints the Hamiltonian terms for the model
	\end{itemize}
	 */
	while ((opt = getopt(argc, argv, "f:s:p:e:o:HBV")) != -1) {
		switch (opt) {
		case 'f':
			filename = optarg;
			break;
		case 's':
			operator_options.site = atoi(optarg);
			break;
		case 'p':
			cmdline_options.precision = atoi(optarg);
			std::cout.precision(cmdline_options.precision);
			std::cerr.precision(cmdline_options.precision);
			break;
		case 'e':
			operator_options.introspect
			    = OptionsForIntrospect::IntrospectEnum::EXPRESSION;
			operator_options.opexpr = optarg;
			break;
		case 'o':
			cmdline_options.solver_options += optarg;
			break;
		case 'B':
			operator_options.introspect
			    = OptionsForIntrospect::IntrospectEnum::MODEL_BASIS;
			break;
		case 'H':
			operator_options.introspect
			    = OptionsForIntrospect::IntrospectEnum::MODEL_HAMILTONIAN;
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

	if (optind < argc) {
		std::cerr << "WARNING: Garbage at end of command line will be ignored\n";
	}

	PsimagLite::String data;
	PsimagLite::InputNg<InputCheck>::Writeable::readFile(data, filename);
	DmrgRunner<double> dmrg_runner(application, data, cmdline_options);

	dmrg_runner.doOneRun(operator_options);
}
