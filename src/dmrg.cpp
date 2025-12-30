#include "BlockDiagonalMatrix.h"
#include "Concurrency.h"
#include "DmrgDriver.h"
#include "Io/IoNg.h"
#include "Provenance.h"
#include "RedirectOutput.hh"
#include "RegisterSignals.h"
#include "RunFinished.h"

typedef PsimagLite::Vector<PsimagLite::String>::Type VectorStringType;

using namespace Dmrg;

typedef PsimagLite::Concurrency ConcurrencyType;

void printLicense(const PsimagLite::PsiApp& app, const OperatorOptions& options)
{
	if (!ConcurrencyType::root() || options.enabled)
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

template <typename MatrixVectorType, typename VectorWithOffsetType>
void mainLoop3(typename MatrixVectorType::ModelType::SuperGeometryType& geometry,
               const ParametersDmrgSolverType&                          dmrgSolverParams,
               InputNgType::Readable&                                   io,
               const OperatorOptions&                                   opOptions)
{
	typedef PsimagLite::ParametersForSolver<typename MatrixVectorType::RealType>
	    ParametersForSolverType;
#ifndef MIN_COMPILE
	if (dmrgSolverParams.options.isSet("ChebyshevSolver")) {
		typedef PsimagLite::ChebyshevSolver<ParametersForSolverType,
		                                    MatrixVectorType,
		                                    typename MatrixVectorType::VectorType>
		    SolverType;
		mainLoop4<SolverType, VectorWithOffsetType>(
		    geometry, dmrgSolverParams, io, opOptions);
	} else {
#else
	{
#endif
		typedef PsimagLite::LanczosSolver<ParametersForSolverType,
		                                  MatrixVectorType,
		                                  typename MatrixVectorType::VectorType>
		    SolverType;
		mainLoop4<SolverType, VectorWithOffsetType>(
		    geometry, dmrgSolverParams, io, opOptions);
	}
}

template <typename MatrixVectorType>
void mainLoop2(InputNgType::Readable&          io,
               const ParametersDmrgSolverType& dmrgSolverParams,
               const OperatorOptions&          opOptions)
{

	typedef typename MatrixVectorType::ModelType        ModelBaseType;
	typedef typename ModelBaseType::QnType              QnType;
	typedef typename ModelBaseType::ModelHelperType     ModelHelperType;
	typedef typename ModelHelperType::SparseElementType SparseElementType;
	typedef SuperGeometry<SparseElementType, InputNgType::Readable, ProgramGlobals>
	    SuperGeometryType;

	SuperGeometryType geometry(io);
	if (dmrgSolverParams.options.isSet("printgeometry"))
		std::cout << geometry;
#ifndef MIN_COMPILE

	if (dmrgSolverParams.options.isSet("vectorwithoffsets")) {
		typedef VectorWithOffsets<SparseElementType, QnType> VectorWithOffsetType;
		mainLoop3<MatrixVectorType, VectorWithOffsetType>(
		    geometry, dmrgSolverParams, io, opOptions);
	} else {
#else
	{
#endif
		typedef VectorWithOffset<SparseElementType, QnType> VectorWithOffsetType;
		mainLoop3<MatrixVectorType, VectorWithOffsetType>(
		    geometry, dmrgSolverParams, io, opOptions);
	}
}

template <typename ComplexOrRealType>
void mainLoop1(InputNgType::Readable&          io,
               const ParametersDmrgSolverType& dmrgSolverParams,
               const OperatorOptions&          opOptions)
{
	typedef SuperGeometry<ComplexOrRealType, InputNgType::Readable, ProgramGlobals>
	                                                          SuperGeometryType;
	typedef PsimagLite::CrsMatrix<ComplexOrRealType>          MySparseMatrix;
	typedef Basis<MySparseMatrix>                             BasisType;
	typedef BasisWithOperators<BasisType>                     BasisWithOperatorsType;
	typedef LeftRightSuper<BasisWithOperatorsType, BasisType> LeftRightSuperType;
	typedef ModelHelperLocal<LeftRightSuperType>              ModelHelperType;
	typedef ModelBase<ModelHelperType,
	                  ParametersDmrgSolverType,
	                  InputNgType::Readable,
	                  SuperGeometryType>
	    ModelBaseType;

#ifndef MIN_COMPILE
	if (dmrgSolverParams.options.isSet("MatrixVectorStored"))
		mainLoop2<MatrixVectorStored<ModelBaseType>>(io, dmrgSolverParams, opOptions);
	else if (dmrgSolverParams.options.isSet("MatrixVectorOnTheFly"))

		mainLoop2<MatrixVectorOnTheFly<ModelBaseType>>(io, dmrgSolverParams, opOptions);
	else
#endif
		mainLoop2<MatrixVectorKron<ModelBaseType>>(io, dmrgSolverParams, opOptions);
}

int main(int argc, char** argv)
{
	PsimagLite::PsiApp application("DMRG++", &argc, &argv, 1);
	InputCheck         inputCheck;
	PsimagLite::String filename = "";
	int                opt      = 0;
	OperatorOptions    options;
	PsimagLite::String strUsage(application.name());
	if (PsimagLite::basename(argv[0]) == "operator")
		options.enabled = true;
	strUsage += " -f filename [-k] [-p precision] [-o solverOptions] [-V] [whatToMeasure]";
	PsimagLite::String sOptions("");
	int                precision    = 0;
	bool               unbuffered   = false;
	SizeType           threadsInCmd = 0;
	bool               versionOnly  = false;
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
	/* PSIDOC OperatorDriver
	 The arguments to the \verb!operator! executable are as follows.
	\begin{itemize}
	 \item[-f] [Mandatory, String] Input to use. The Model= line is
	very important in input.inp.

	\item[-e] [Mandatory unless -l, String] OperatorExpression; see manual

	\item[-s] [Optional, Integer] \emph{Deprecated. Use -e.}
	Site for operator.
	Meaningful only for Models where
	the Hilbert space depends on the site (different kinds of atoms).
	Defaults to 0.

	\item[-l] [Mandatory unless -e, String] \emph{Deprecated. Use -e.}
	The label or name for this operator.
	This is model dependent. For example to obtain $c_{\uparrow}$ for
	the Hubbard model, say \begin{verbatim}
	./operator -l c -f input.inp\end{verbatim}
	See the function naturalOperator for each Model.

	\item[-B] [Optional] Prints the basis and all operators for the model

	\item[-H] [Optional] Prints the Hamiltonian terms for the model

	\item[-d] [Optional, Integer] \emph{Deprecated. Use -e.}
	Degree of freedom (spin, orbital or
	combination of both) to use. This is model dependent. For example to
	obtain $c_\downarrow$ for the Hubbard model, say
	\begin{verbatim}./operator -l c -d 1 -f input.inp\end{verbatim}
	See the function naturalOperator for each Model. Defaults to 0.

	\item[-t] [Optional, Void] \emph{Note: When using -e, it transposes
	the whole expression.}
	Transpose the operator. For example to
	obtain $c^\dagger_\uparrow$ for a Hubbard model, say
	\begin{verbatim}./operator -l c -t -f input.inp\end{verbatim}
	\end{itemize}
	 */
	while ((opt = getopt(argc, argv, "f:s:l:d:p:e:o:S:tkBHUV")) != -1) {
		switch (opt) {
		case 'f':
			filename = optarg;
			break;
		case 's':
			options.site = atoi(optarg);
			break;
		case 'l':
			options.label = optarg;
			break;
		case 'd':
			options.dof = atoi(optarg);
			break;
		case 't':
			options.transpose = true;
			break;
		case 'p':
			precision = atoi(optarg);
			std::cout.precision(precision);
			std::cerr.precision(precision);
			break;
		case 'e':
			options.hasOperatorExpression = true;
			options.opexpr                = optarg;
			break;
		case 'o':
			sOptions += optarg;
			break;
		case 'S':
			threadsInCmd = atoi(optarg);
			break;
		case 'B':
			options.label = "B";
			break;
		case 'H':
			options.label = "H";
		case 'U':
			unbuffered = true;
			break;
		case 'V':
			versionOnly   = true;
			options.label = "-";
			break;
		default:
			inputCheck.usageMain(strUsage);
			return 1;
		}
	}

	// sanity checks here
	if (filename == "" && !versionOnly) {
		inputCheck.usageMain(strUsage);
		return 1;
	}

	PsimagLite::String insitu = (optind < argc) ? argv[optind] : "";

	if (!options.enabled && options.label != "-") {
		bool queryOnly = (options.label == "?");
		if (options.label == "" || options.label == "?") {
			options.label = ProgramGlobals::coutName(filename);
			if (queryOnly) {
				std::cout << options.label << "\n";
				return 0;
			}
		}
	}

	// print license
	if (versionOnly) {
		printLicense(application, options);
		return 0;
	}

	InputNgType::Writeable ioWriteable(filename, inputCheck);
	InputNgType::Readable  io(ioWriteable);

	ParametersDmrgSolverType dmrgSolverParams(io, sOptions, false);

	if (dmrgSolverParams.options.isSet("hd5DontPrint"))
		PsimagLite::IoNg::dontPrintDebug();

	if (dmrgSolverParams.options.isSet("addPidToOutputName"))
		options.label += "." + ttos(getpid());

	if (threadsInCmd > 0)
		dmrgSolverParams.nthreads = threadsInCmd;
	if (precision > 0)
		dmrgSolverParams.precision = precision;

	bool echoInput = false;
	if (!options.enabled && options.label != "-") {

		RunFinished runFinished(dmrgSolverParams.options.isSet("noClobber"));
		if (runFinished.OK(options.label.c_str())) {
			runFinished.printTermination(options.label.c_str());
			return 1;
		}

		echoInput = true;

		PsimagLite::RedirectOutput::setAppName(application.name(),
		                                       Provenance::logo(application.name()));

		std::ios_base::openmode open_mode
		    = (dmrgSolverParams.autoRestart) ? std::ofstream::app : std::ofstream::out;

		PsimagLite::RedirectOutput::doIt(options.label, open_mode, unbuffered);
	}

	if (dmrgSolverParams.autoRestart) {
		std::cout << "\nAutoRestart possible\n";
		std::cerr << "AutoRestart possible\n";
	}

	printLicense(application, options);

	application.printCmdLine(std::cout);
	if (echoInput)
		application.echoBase64(std::cout, filename);

	if (insitu != "")
		dmrgSolverParams.insitu = insitu;
	if (dmrgSolverParams.options.isSet("minimizeDisk"))
		dmrgSolverParams.options += ",noSaveWft,noSaveStacks,noSaveData";

	bool setAffinities = (dmrgSolverParams.options.isSet("setAffinities"));

	SizeType threadsStackSize = 0;
	try {
		io.readline(threadsStackSize, "ThreadsStackSize=");
	} catch (std::exception&) { }

	PsimagLite::CodeSectionParams codeSection(
	    dmrgSolverParams.nthreads, dmrgSolverParams.nthreads2, setAffinities, threadsStackSize);
	ConcurrencyType::setOptions(codeSection);

	registerSignals();

	bool isComplex = (dmrgSolverParams.options.isSet("useComplex")
	                  || dmrgSolverParams.options.isSet("TimeStepTargeting"));

	if (isComplex)
		mainLoop1<std::complex<RealType>>(io, dmrgSolverParams, options);
	else
		mainLoop1<RealType>(io, dmrgSolverParams, options);
}
