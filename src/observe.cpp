#include "ObserveDriver.h"
#include "RedirectOutput.hh"

using namespace Dmrg;

template <typename T>
bool atLeastOneLoopWithBit0Set(const T& fl)
{
	const SizeType n = fl.size();
	for (SizeType i = 0; i < n; ++i)
		if (fl[i].wants("save"))
			return true;

	return false;
}
template <typename GeometryType,
    typename ModelHelperType,
    typename VectorWithOffsetType>
void mainLoop(GeometryType& geometry,
    InputNgType::Readable& io,
    const ParametersDmrgSolverType& params,
    const std::string& list)
{
	typedef typename VectorWithOffsetType::value_type ComplexOrRealType;

	typedef ModelBase<ModelHelperType,
	    ParametersDmrgSolverType,
	    InputNgType::Readable,
	    GeometryType>
	    ModelBaseType;

	SizeType orbitals = 1.0;
	try {
		io.readline(orbitals, "Orbitals=");
	} catch (std::exception&) {
	}

	ModelSelector<ModelBaseType> modelSelector(params.model);
	const ModelBaseType& model = modelSelector(params, io, geometry);

	const std::string& datafile = params.filename;
	IoInputType dataIo(datafile);

	bool iscomplex = false;
	dataIo.read(iscomplex, "IsComplex");

	if (iscomplex != PsimagLite::IsComplexNumber<ComplexOrRealType>::True)
		err("Previous run was complex and this one is not (or viceversa)\n");

	while (!observeOneFullSweep<VectorWithOffsetType, ModelBaseType>(dataIo, model, list, orbitals))
		;
}

template <typename GeometryType,
    template <typename> class ModelHelperTemplate,
    typename MySparseMatrix>
void mainLoop1(GeometryType& geometry,
    InputNgType::Readable& io,
    const ParametersDmrgSolverType& params,
    const std::string& list)
{
	typedef Basis<MySparseMatrix> BasisType;
	typedef BasisWithOperators<BasisType> BasisWithOperatorsType;
	typedef LeftRightSuper<BasisWithOperatorsType, BasisType> LeftRightSuperType;
	typedef ModelHelperTemplate<LeftRightSuperType> ModelHelperType;
	typedef typename MySparseMatrix::value_type ComplexOrRealType;
	typedef Qn QnType;

	if (params.options.isSet("vectorwithoffsets")) {
		typedef VectorWithOffsets<ComplexOrRealType, QnType> VectorWithOffsetType;
		mainLoop<GeometryType, ModelHelperType, VectorWithOffsetType>(geometry, io, params, list);
	} else {
		typedef VectorWithOffset<ComplexOrRealType, QnType> VectorWithOffsetType;
		mainLoop<GeometryType, ModelHelperType, VectorWithOffsetType>(geometry, io, params, list);
	}
}

template <typename MySparseMatrix>
void mainLoop0(InputNgType::Readable& io,
    ParametersDmrgSolverType& dmrgSolverParams,
    InputCheck& inputCheck,
    const std::string& list)
{
	typedef typename MySparseMatrix::value_type ComplexOrRealType;
	typedef Dmrg::SuperGeometry<ComplexOrRealType,
	    InputNgType::Readable,
	    ProgramGlobals>
	    SuperGeometryType;

	SuperGeometryType superGeometry(io);
	int tmp = 0;
	try {
		io.readline(tmp, "UseSu2Symmetry=");
	} catch (std::exception&) {
	}

	bool su2 = (tmp > 0);

	if (su2) {
		err("SU(2) no longer supported\n");
	} else {
		mainLoop1<SuperGeometryType, ModelHelperLocal, MySparseMatrix>(superGeometry, io, dmrgSolverParams, list);
	}
}

void usage(const std::string& name)
{
	std::cerr << "USAGE is " << name << " -f filename [-l output_file] [-p precision] ";
	std::cerr << " [-F fileoption] [-U] [-V] whatToMeasure\n";
}

/* PSIDOC ObserveDriver
 The observe driver may be used to measure after the main dmrg run has finished.
The observe driver is used as follows,
\begin{verbatim}
./observe -f input.inp whatToMeasure
\end{verbatim}
  The command line arguments of observe are the following.
\begin{itemize}
\item[-f] [Mandatory, String] Input to use. Files
 referred to by \verb!OutputFile=! are now inputs, and
must be present.
\item[whatToMeasure] {[}Mandatory, String{]} What to measure post process.
This is a comma-separated list of braket specifications.
Braket specifications can be bare or dressed, and are explained elsewhere.
\item[-l] {[}Optional, String{]} Redirect std::cout to this file.
\item[-o] {[}Optional, String{]} Extra options for SolverOptions
\item[-p] [Optional, Integer] Digits of precision for printing.
\item[-F] [Optional, string] TBW
\item[-S] [Optional, number] Ignore the Threads= line if present in the input, and run with
Threads=number
\item[-U] [Optional] Make cout output unbuffered
\item[-V] [Optional] Print version and exit
\end{itemize}
  */
int main(int argc, char** argv)
{
	using namespace Dmrg;
	PsimagLite::PsiApp application("observe", &argc, &argv, 1);
	std::string filename;
	std::string filesOption;
	std::string output_filename;
	std::string sOptions(",observe");
	int opt = 0;
	int precision = 0;
	SizeType threadsInCmd = 0;
	bool versionOnly = false;
	bool unbuffered = false;

	while ((opt = getopt(argc, argv, "f:l:o:p:F:S:UV")) != -1) {
		switch (opt) {
		case 'f':
			filename = optarg;
			break;
		case 'l':
			output_filename = optarg;
			break;
		case 'p':
			precision = atoi(optarg);
			std::cout.precision(precision);
			std::cerr.precision(precision);
			break;
		case 'o':
			sOptions += optarg;
			break;
		case 'F':
			filesOption = optarg;
			break;
		case 'S':
			threadsInCmd = atoi(optarg);
			break;
		case 'U':
			unbuffered = true;
			break;
		case 'V':
			versionOnly = true;
			break;
		default:
			usage(application.name());
			return 1;
		}
	}

	std::string list = (optind < argc) ? argv[optind] : "";

	// sanity checks here
	if (filename == "" || (filesOption != "keep" && filesOption != "")) {
		if (!versionOnly) {
			usage(application.name());
			return 1;
		}
	}

	if (!output_filename.empty()) {
		bool queryOnly = (output_filename == "?");
		if (queryOnly) {
			output_filename = ProgramGlobals::coutName(filename, application.name());
			std::cout << output_filename << "\n";
			return 0;
		}

		bool redirect = (output_filename == "+r");
		if (redirect) {
			output_filename = ProgramGlobals::coutName(filename, application.name());
		}

		PsimagLite::RedirectOutput::setAppName(application.name(),
		    Provenance::logo(application.name()));

		std::ios_base::openmode open_mode = std::ofstream::out;

		PsimagLite::RedirectOutput::doIt(output_filename, open_mode, unbuffered);
		application.echoBase64(std::cout, filename);
	}

	typedef PsimagLite::Concurrency ConcurrencyType;

	// print license
	if (ConcurrencyType::root()) {
		Provenance provenance;
		std::cout << provenance;
		std::cout << Provenance::logo(application.name()) << "\n";
		application.checkMicroArch(std::cout, Provenance::compiledMicroArch());
	}

	if (versionOnly)
		return 0;

	application.printCmdLine(std::cout);

	// Setup the Geometry
	InputCheck inputCheck;
	InputFromDataOrNot<InputCheck> inputFromDataOrNot(filename, inputCheck, false);
	InputNgType::Readable io(inputFromDataOrNot.ioWriteable());

	ParametersDmrgSolverType dmrgSolverParams(io, sOptions, false, true);

	if (dmrgSolverParams.options.isSet("hd5DontPrint"))
		PsimagLite::IoNg::dontPrintDebug();

	if (threadsInCmd > 0)
		dmrgSolverParams.nthreads = threadsInCmd;
	if (precision > 0)
		dmrgSolverParams.precision = precision;

	bool setAffinities = (dmrgSolverParams.options.isSet("setAffinities"));

	SizeType threadsStackSize = 0;
	try {
		io.readline(threadsStackSize, "ThreadsStackSize=");
	} catch (std::exception&) {
	}

	PsimagLite::CodeSectionParams codeSectionParams(dmrgSolverParams.nthreads,
	    dmrgSolverParams.nthreads2,
	    setAffinities,
	    threadsStackSize);
	ConcurrencyType::setOptions(codeSectionParams);

	if (!atLeastOneLoopWithBit0Set(dmrgSolverParams.finiteLoop))
		err("FATAL: At least one loop must have bit 0 set for observe to work\n");

	bool isComplex = (dmrgSolverParams.options.isSet("useComplex") || dmrgSolverParams.options.isSet("TimeStepTargeting"));

	if (isComplex) {
		mainLoop0<MySparseMatrixComplex>(io, dmrgSolverParams, inputCheck, list);
	} else {
		mainLoop0<MySparseMatrixReal>(io, dmrgSolverParams, inputCheck, list);
	}

} // main
