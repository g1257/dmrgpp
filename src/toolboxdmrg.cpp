#include "InputCheck.h"
#include "InputFromDataOrNot.h"
#include "InputNg.h"
#include "ParametersDmrgSolver.h"
#include "ProgramGlobals.h"
#include "PsimagLite.h"
#include "Qn.h"
#include "ToolBox.h"
#include <iostream>

#ifndef USE_FLOAT
using RealType = double;
#else
using RealType = float;
#endif
using InputNgType = PsimagLite::InputNg<Dmrg::InputCheck>;
using ParametersDmrgSolverType
    = Dmrg::ParametersDmrgSolver<RealType, InputNgType::Readable, Dmrg::Qn>;
using ConcurrencyType = PsimagLite::Concurrency;

void usage(const PsimagLite::String& name)
{
	std::cerr << "USAGE is " << name << " -f filename -a action [-s] [-p precision]\n";
}

struct ToolOptions {

	ToolOptions()
	    : extraOptions("lowest eigenvalue")
	    , shortoption(false)
	{ }

	PsimagLite::String filename;
	PsimagLite::String action;
	PsimagLite::String extraOptions;
	bool               shortoption;
};

template <typename ComplexOrRealType>
void main1(InputNgType::Readable&          io,
           PsimagLite::PsiApp              application,
           const ParametersDmrgSolverType& dmrgSolverParams,
           const ToolOptions&              toolOptions)
{
	using GeometryType
	    = PsimagLite::Geometry<ComplexOrRealType, InputNgType::Readable, Dmrg::ProgramGlobals>;
	GeometryType geometry(io);

	using ToolBoxType = Dmrg::ToolBox<ParametersDmrgSolverType, GeometryType>;
	ConcurrencyType::codeSectionParams.npthreads = dmrgSolverParams.nthreads;
	PsimagLite::String label
	    = (toolOptions.action == "energies") ? "lowest" : toolOptions.extraOptions;
	typename ToolBoxType::ParametersForGrepType params(label, toolOptions.shortoption);
	typename ToolBoxType::ActionEnum act = ToolBoxType::actionCanonical(toolOptions.action);
	if (act == ToolBoxType::ACTION_GREP) {
		ToolBoxType::printGrep(toolOptions.filename, params);
	} else if (act == ToolBoxType::ACTION_INPUT) {
		std::cout << io.data() << "\n";
	} else if (act == ToolBoxType::ACTION_ANALYSIS) {
		PsimagLite::String str("Analyzing ");
		str += toolOptions.filename;
		std::cout << str << "\n";
		ToolBoxType::analize(dmrgSolverParams, geometry, toolOptions.extraOptions);
	} else {
		std::cerr << application.name();
		std::cerr << ": Unknown action " << toolOptions.action << "\n";
		std::cerr << "\tSupported actions are " << ToolBoxType::actions() << "\n";
	}
}

/* PSIDOC ToolboxDriver
 The command line arguments of toolboxdmrg are the following.
  \begin{itemize}
  \item[-f] {[}Mandatory, String{]} Input to use. Files
 referred to by \verb!OutputFile=! are now inputs, and
must be present.
  \item[-a] {[}Mandatory, String{]} Action, see below.
  \item[-E] {[}Optional, String{]} Extra options, see below.
  \item[-o] {[}Optional, String{]} Extra options for SolverOptions
  \item[-s] {[}Optional, no argument needed{]} Short option.
  \item[-p] [Optional, Integer] Digits of precision for printing.
 \item[-V] [Optional] Print version and exit
  \end{itemize}
*/
int main(int argc, char** argv)
{
	using namespace Dmrg;
	PsimagLite::PsiApp application("toolboxdmrg", &argc, &argv, 1);
	ToolOptions        toolOptions;
	int                opt         = 0;
	int                precision   = 0;
	bool               versionOnly = false;
	PsimagLite::String sOptions;
	while ((opt = getopt(argc, argv, "f:p:a:E:o:sV")) != -1) {
		switch (opt) {
		case 'f':
			toolOptions.filename = optarg;
			break;
		case 'p':
			precision = atoi(optarg);
			std::cout.precision(precision);
			std::cerr.precision(precision);
			break;
		case 'a':
			toolOptions.action = optarg;
			break;
		case 'E':
			toolOptions.extraOptions = optarg;
			break;
		case 'o':
			sOptions += optarg;
			break;
		case 's':
			toolOptions.shortoption = true;
			break;
		case 'V':
			versionOnly = true;
			break;
		default:
			usage(application.name());
			return 1;
		}
	}

	// sanity checks here
	if (toolOptions.filename == "" || toolOptions.action == "") {
		if (!versionOnly) {
			usage(application.name());
			return 1;
		}
	}

	// print license
	if (ConcurrencyType::root()) {
		std::cerr << ProgramGlobals::license;
		Provenance provenance;
		std::cout << provenance;
	}

	if (versionOnly)
		return 0;

	InputCheck inputCheck;

	bool               filenameIsCout = false;
	PsimagLite::String dotCout        = ".cout";
	size_t             pos            = toolOptions.filename.find(dotCout);
	if (pos != PsimagLite::String::npos
	    && pos + dotCout.size() == toolOptions.filename.size()) {
		filenameIsCout = true;
	}

	if (toolOptions.action == "files") {
		if (toolOptions.extraOptions == "")
			toolOptions.extraOptions = "list";
		inputCheck.checkFileOptions(toolOptions.extraOptions);
	}

	InputFromDataOrNot<InputCheck> inputFromDataOrNot(
	    toolOptions.filename, inputCheck, filenameIsCout);
	InputNgType::Readable io(inputFromDataOrNot.ioWriteable());

	//! Read the parameters for this run
	bool                     earlyExit = true;
	ParametersDmrgSolverType dmrgSolverParams(io, sOptions, earlyExit);

	if (precision > 0)
		dmrgSolverParams.precision = precision;

	if (dmrgSolverParams.options.isSet("useComplex"))
		main1<std::complex<RealType>>(io, application, dmrgSolverParams, toolOptions);
	else
		main1<RealType>(io, application, dmrgSolverParams, toolOptions);
} // main
