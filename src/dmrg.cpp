#include "BlockDiagonalMatrix.h"
#include "IoSimple.h"
#include "Concurrency.h"
#include "Provenance.h"
#include "RegisterSignals.h"
#include "ArchiveFiles.h"
#include "DmrgDriver.h"

typedef PsimagLite::Vector<PsimagLite::String>::Type VectorStringType;
typedef  PsimagLite::CrsMatrix<std::complex<RealType> > MySparseMatrixComplex;
typedef  PsimagLite::CrsMatrix<RealType> MySparseMatrixReal;

using namespace Dmrg;

typedef ArchiveFiles<ParametersDmrgSolverType> ArchiveFilesType;

std::streambuf *GlobalCoutBuffer = 0;
std::ofstream GlobalCoutStream;

void usageOperator()
{
	std::cerr<<"USAGE is operator -f filename ";
	std::cerr<<"-l label [-d dof] [-s site] [-t]\n";
}

void restoreCoutBuffer()
{
	if (GlobalCoutBuffer == 0) return;
	GlobalCoutStream.close();
	std::cout.rdbuf(GlobalCoutBuffer);
}

template<typename MatrixVectorType, typename VectorWithOffsetType>
void mainLoop3(typename MatrixVectorType::ModelType::GeometryType& geometry,
               const ParametersDmrgSolverType& dmrgSolverParams,
               InputNgType::Readable& io,
               const OperatorOptions& opOptions,
               PsimagLite::String targeting)
{
	typedef PsimagLite::ParametersForSolver<typename MatrixVectorType::RealType>
	        ParametersForSolverType;
	if (dmrgSolverParams.options.find("ChebyshevSolver")!=PsimagLite::String::npos) {
		typedef PsimagLite::ChebyshevSolver<ParametersForSolverType,
		        MatrixVectorType, typename MatrixVectorType::VectorType> SolverType;
		mainLoop4<SolverType,VectorWithOffsetType>(geometry,
		                                           dmrgSolverParams,
		                                           io,
		                                           opOptions,
		                                           targeting);
	} else {
		typedef PsimagLite::LanczosSolver<ParametersForSolverType,
		        MatrixVectorType, typename MatrixVectorType::VectorType> SolverType;
		mainLoop4<SolverType,VectorWithOffsetType>(geometry,
		                                           dmrgSolverParams,
		                                           io,
		                                           opOptions,
		                                           targeting);
	}
}

template<typename MatrixVectorType>
void mainLoop2(typename MatrixVectorType::ModelType::GeometryType& geometry,
               const ParametersDmrgSolverType& dmrgSolverParams,
               InputNgType::Readable& io,
               const OperatorOptions& opOptions,
               PsimagLite::String targeting)
{
	typedef typename MatrixVectorType::ComplexOrRealType ComplexOrRealType;

	if (dmrgSolverParams.options.find("vectorwithoffsets")!=PsimagLite::String::npos) {
		typedef VectorWithOffsets<ComplexOrRealType> VectorWithOffsetType;
		mainLoop3<MatrixVectorType,VectorWithOffsetType>(geometry,
		                                                 dmrgSolverParams,
		                                                 io,
		                                                 opOptions,
		                                                 targeting);
	} else {
		typedef VectorWithOffset<ComplexOrRealType> VectorWithOffsetType;
		mainLoop3<MatrixVectorType,VectorWithOffsetType>(geometry,
		                                                 dmrgSolverParams,
		                                                 io,
		                                                 opOptions,
		                                                 targeting);
	}
}

template<typename GeometryType,
         template<typename> class ModelHelperTemplate,
         typename MySparseMatrix,
         typename CvectorSizeType>
void mainLoop1(GeometryType& geometry,
               const ParametersDmrgSolverType& dmrgSolverParams,
               InputNgType::Readable& io,
               const OperatorOptions& opOptions,
               PsimagLite::String targeting)
{
	typedef Basis<MySparseMatrix, CvectorSizeType> BasisType;
	typedef Operators<BasisType> OperatorsType;
	typedef BasisWithOperators<OperatorsType> BasisWithOperatorsType;
	typedef LeftRightSuper<BasisWithOperatorsType,BasisType> LeftRightSuperType;
	typedef ModelHelperTemplate<LeftRightSuperType> ModelHelperType;
	typedef ModelBase<ModelHelperType,
	        ParametersDmrgSolverType,
	        InputNgType::Readable,
	        GeometryType> ModelBaseType;

	if (dmrgSolverParams.options.find("MatrixVectorStored")!=PsimagLite::String::npos) {
		mainLoop2<MatrixVectorStored<ModelBaseType> >(geometry,
		                                              dmrgSolverParams,
		                                              io,
		                                              opOptions,
		                                              targeting);
	} else if (dmrgSolverParams.options.find("MatrixVectorKron")!=PsimagLite::String::npos) {
		if (ModelHelperType::isSu2())
			throw PsimagLite::RuntimeError("MatrixVectorKron does not support SU(2)\n");

		mainLoop2<MatrixVectorKron<ModelBaseType> >(geometry,
		                                            dmrgSolverParams,
		                                            io,
		                                            opOptions,
		                                            targeting);
	} else {
		mainLoop2<MatrixVectorOnTheFly<ModelBaseType> >(geometry,
		                                                dmrgSolverParams,
		                                                io,
		                                                opOptions,
		                                                targeting);
	}
}

template<typename MySparseMatrix, typename CvectorSizeType>
void mainLoop0(InputNgType::Readable& io,
               const ParametersDmrgSolverType& dmrgSolverParams,
               PsimagLite::String targeting,
               const OperatorOptions& opOptions)
{
	typedef typename MySparseMatrix::value_type ComplexOrRealType;
	typedef PsimagLite::Geometry<ComplexOrRealType,
	        InputNgType::Readable,
	        ProgramGlobals> GeometryType;

	GeometryType geometry(io);
	if (dmrgSolverParams.options.find("printgeometry") != PsimagLite::String::npos)
		std::cout<<geometry;

	int tmp = 0;
	try {
		io.readline(tmp,"UseSu2Symmetry=");
	} catch (std::exception&) {}

	bool su2 = (tmp > 0);

	if (targeting != "GroundStateTargetting" && su2) {
		PsimagLite::String str("SU(2) supports only GroundStateTargetting (sorry!)\n");
		throw PsimagLite::RuntimeError(str);
	}

	if (su2) {
		mainLoop1<GeometryType,ModelHelperSu2,MySparseMatrix,CvectorSizeType>(geometry,
		                                                                      dmrgSolverParams,
		                                                                      io,
		                                                                      opOptions,
		                                                                      targeting);
		return;
	}

	if (dmrgSolverParams.options.find("useComplex") != PsimagLite::String::npos &&
	        targeting != "TimeStepTargetting" &&
	        targeting != "GroundStateTargetting" &&
		targeting != "TargetingCorrelations" &&
		targeting != "CorrectionTargetting") {
		PsimagLite::String str("SolverOptions=useComplex not allowed for ");
		str += targeting + "\n";
		throw PsimagLite::RuntimeError(str);
	}

	mainLoop1<GeometryType,ModelHelperLocal,MySparseMatrix,CvectorSizeType>(geometry,
	                                                                        dmrgSolverParams,
	                                                                        io,
	                                                                        opOptions,
	                                                                        targeting);
}

int main(int argc, char *argv[])
{
	PsimagLite::PsiApp application("DMRG++",&argc,&argv,1);
	typedef PsimagLite::Concurrency ConcurrencyType;
	InputCheck inputCheck;
	PsimagLite::String filename="";
	int opt = 0;
	OperatorOptions options;
	PsimagLite::String strUsage(application.name());
	if (utils::basename(argv[0]) == "operator") options.enabled = true;
	strUsage += " -f filename [-k] [-p precision] [-V] [whatToMeasure]";
	int precision = 6;
	bool keepFiles = false;
	bool versionOnly = false;
	/* PSIDOC DmrgDriver
	  \begin{itemize}
	  \item[-f] {[}Mandatory, String{]} Input to use.
	  \item[-p] [Optional, Integer] Digits of precision for printing.
	  \item[whatToMeasure] {[}Optional, String{]} What to measure in-situ
	  \item[-l] {[}Optional, String{]} Without this option std::cout is redirected
	  to a file.
	  This option with the string ``?'' prints name of such log file.
	  This option with the string ``-'' writes std::cout to terminal.
	  In other cases, string is the name of the file to redirect std::cout to.
	 \item[-k] [Optional] Keep untar files
	  \end{itemize}
	 */
	/* PSIDOC OperatorDriver
	 The arguments to the \verb!operator! executable are as follows.
	\begin{itemize}
	 \item[-f] [Mandatory, String] Input to use. The Model= line is
	very important in input.inp.

	\item[-s] [Optional, Integer] Site for operator.
	Meaningful only for Models where
	the Hilbert space depends on the site (different kinds of atoms).
	Defaults to 0.

	\item[-l] [Mandatory, String] The label or name for this operator.
	This is model dependent. For example to obtain $c_{\uparrow}$ for
	the Hubbard model, say \begin{verbatim}
	./operator -l c -f input.inp\end{verbatim}
	See the function naturalOperator for each Model.

	\item[-d] [Optional, Integer] Degree of freedom (spin, orbital or
	combination of both) to use. This is model dependent. For example to
	obtain $c_\downarrow$ for the Hubbard model, say
	\begin{verbatim}./operator -l c -d 1 -f input.inp\end{verbatim}
	See the function naturalOperator for each Model. Defaults to 0.

	\item[-t] [Optional, Void] Transpose the operator. For example to
	obtain $c^\dagger_\uparrow$ for a Hubbard model, say
	\begin{verbatim}./operator -l c -t -f input.inp\end{verbatim}
	\end{itemize}
	 */
	while ((opt = getopt(argc, argv,"f:s:l:d:p:tkV")) != -1) {
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
		case 'k':
			keepFiles = true;
			break;
		case 'p':
			precision = atoi(optarg);
			std::cout.precision(precision);
			std::cerr.precision(precision);
			break;
		case 'V':
			versionOnly = true;
			options.label = "-";
			break;
		default:
			inputCheck.usageMain(strUsage);
			return 1;
		}
	}

	// sanity checks here
	if (filename=="" && !versionOnly) {
		inputCheck.usageMain(strUsage);
		return 1;
	}

	PsimagLite::String insitu = (optind < argc) ? argv[optind] : "";

	if (!options.enabled && options.label != "-") {
		bool queryOnly = (options.label == "?");
		if (options.label == "" || options.label == "?") {
			options.label = ArchiveFilesType::coutName(filename);
			if (queryOnly) {
				std::cout<<options.label<<"\n";
				return 0;
			}
		}

		GlobalCoutStream.open(options.label.c_str());
		if (!GlobalCoutStream || GlobalCoutStream.bad()
		        || !GlobalCoutStream.good()) {
			PsimagLite::String str(application.name());
			str += ": Could not redirect std::cout to " + options.label + "\n";
			throw PsimagLite::RuntimeError(str);
		}

		std::cerr<<application.name()<<"\x1b[38;5;124m";
		std::cerr<<" [features] "<<PsimagLite::AnsiColor::reset;
		std::cerr<<"Standard output sent to ";
		std::cerr<<options.label<<"\n";
		std::cerr.flush();
		GlobalCoutBuffer = std::cout.rdbuf(); //save old buf
		std::cout.rdbuf(GlobalCoutStream.rdbuf()); //redirect std::cout to file
		atexit(restoreCoutBuffer);
	}

	// print license
	if (ConcurrencyType::root() && !options.enabled) {
		std::cout<<ProgramGlobals::license;
		Provenance provenance;
		std::cout<<provenance;
	}

	if (versionOnly) return 0;

	InputNgType::Writeable ioWriteable(filename,inputCheck);
	InputNgType::Readable io(ioWriteable);

	ParametersDmrgSolverType dmrgSolverParams(io, false);

	ArchiveFilesType af(dmrgSolverParams,filename,options.enabled,options.label);

	if (insitu!="") dmrgSolverParams.insitu = insitu;
	if (dmrgSolverParams.options.find("minimizeDisk") != PsimagLite::String::npos)
		dmrgSolverParams.options += ",noSaveWft,noSaveStacks,noSaveData";

	bool setAffinities = (dmrgSolverParams.options.find("setAffinities")
                                  != PsimagLite::String::npos);
	ConcurrencyType::setOptions(dmrgSolverParams.nthreads, setAffinities);

	registerSignals();

	PsimagLite::String targeting = inputCheck.getTargeting(dmrgSolverParams.options);

	bool isComplex = (dmrgSolverParams.options.find("useComplex") != PsimagLite::String::npos);
	if (targeting=="TimeStepTargetting") isComplex = true;
	if (isComplex) {
		mainLoop0<MySparseMatrixComplex, CvectorSizeType>(io,dmrgSolverParams,targeting,options);
	} else {
		mainLoop0<MySparseMatrixReal, CvectorSizeType>(io,dmrgSolverParams,targeting,options);
	}

	if (options.enabled) return 0;

	af.deletePackedFiles();
	if (!keepFiles)
		ArchiveFilesType::staticDelete();
}

