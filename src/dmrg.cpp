#ifndef USE_FLOAT
typedef double MatrixElementType;
#else
typedef float MatrixElementType;
#endif

#include "CrsMatrix.h"
#include "LanczosSolver.h"
#include "ChebyshevSolver.h"
#include "BlockMatrix.h"
#include "DmrgSolver.h"
#include "IoSimple.h"
#include "Operators.h"
#include "Concurrency.h"
#include "Geometry/Geometry.h"
#include "ModelHelperLocal.h"
#include "ModelHelperSu2.h"
#include "MatrixVectorOnTheFly.h"
#include "MatrixVectorStored.h"
#include "MatrixVectorKron.h"
#include "TargetingGroundState.h"
#include "TargetingTimeStep.h"
#include "TargetingDynamic.h"
#include "TargetingAdaptiveDynamic.h"
#include "TargetingCorrection.h"
#include "TargetingCorrectionVector.h"
#include "MettsTargetting.h"
#include "VectorWithOffset.h"
#include "VectorWithOffsets.h"
#include "BasisWithOperators.h"
#include "LeftRightSuper.h"
#include "Provenance.h"
#include "InputNg.h"
#include "InputCheck.h"
#include "ModelSelector.h"
#include "RegisterSignals.h"
#include "ArchiveFiles.h"

typedef PsimagLite::Vector<PsimagLite::String>::Type VectorStringType;
typedef  PsimagLite::CrsMatrix<std::complex<MatrixElementType> > MySparseMatrixComplex;
typedef  PsimagLite::CrsMatrix<MatrixElementType> MySparseMatrixReal;
#ifdef USE_COMPLEX
typedef MySparseMatrixComplex MySparseMatrixC;
#else
typedef MySparseMatrixReal MySparseMatrixC;
#endif

using namespace Dmrg;

typedef PsimagLite::InputNg<InputCheck> InputNgType;
typedef ParametersDmrgSolver<MatrixElementType,InputNgType::Readable> ParametersDmrgSolverType;
typedef ArchiveFiles<ParametersDmrgSolverType> ArchiveFilesType;

struct OperatorOptions {

	OperatorOptions()
	    : site(0),
	      dof(0),
	      label(""),
	      fermionicSign(0),
	      transpose(false),
	      enabled(false)
	{}

	SizeType site;
	SizeType dof;
	PsimagLite::String label;
	int fermionicSign;
	bool transpose;
	bool enabled;
};

std::streambuf *GlobalCoutBuffer = 0;
std::ofstream GlobalCoutStream;

void restoreCoutBuffer()
{
	if (GlobalCoutBuffer == 0) return;
	GlobalCoutStream.close();
	std::cout.rdbuf(GlobalCoutBuffer);
}

void usageOperator()
{
	std::cerr<<"USAGE is operator -f filename -F ";
	std::cerr<<"fermionicSign -l label [-d dof] [-s site] [-t]\n";
}

template<typename ModelBaseType>
void operatorDriver(const ModelBaseType& model, const OperatorOptions& obsOptions)
{
	typedef typename ModelBaseType::ModelHelperType ModelHelperType;
	typedef typename ModelHelperType::OperatorsType OperatorsType;
	typedef typename OperatorsType::OperatorType OperatorType;
	typedef typename ModelHelperType::SparseElementType SparseElementType;
	typedef PsimagLite::Matrix<SparseElementType> MatrixType;
	typedef std::pair<SizeType,SizeType> PairType;

	if (obsOptions.label=="" || obsOptions.fermionicSign == 0) {
		usageOperator();
		return;
	}

	MatrixType opC = model.naturalOperator(obsOptions.label,obsOptions.site,obsOptions.dof);
	std::cerr<<"#label="<<obsOptions.label<<" site="<<obsOptions.site;
	std::cerr<<" dof="<<obsOptions.dof<<"\n";

	Su2Related su2Related;

	MatrixType opC2;
	if (obsOptions.transpose)
		transposeConjugate(opC2,opC);

	OperatorType opC3((obsOptions.transpose) ? opC2 : opC,
	                  obsOptions.fermionicSign,
	                  PairType(0,0),
	                  1,
	                  su2Related);

	opC3.save(std::cout);
}

template<typename GeometryType,
         typename TargettingType>
void mainLoop3(GeometryType& geometry,
               const ParametersDmrgSolverType& dmrgSolverParams,
               InputNgType::Readable& io,
               const OperatorOptions& opOptions)
{
	typedef typename TargettingType::TargetParamsType TargetParamsType;
	typedef typename TargettingType::MatrixVectorType::ModelType ModelBaseType;

	//! Setup the Model
	ModelSelector<ModelBaseType> modelSelector(dmrgSolverParams.model);
	const ModelBaseType& model = modelSelector(dmrgSolverParams,io,geometry);

	if (opOptions.enabled) {
		operatorDriver(model,opOptions);
		return;
	}

	//! Read TimeEvolution if applicable:
	TargetParamsType tsp(io,model);

	//! Setup the dmrg solver:
	typedef DmrgSolver<TargettingType> SolverType;
	SolverType dmrgSolver(model,tsp,io);

	//! Calculate observables:
	dmrgSolver.main(geometry);
}

template<typename GeometryType,
         typename MatrixVectorType,
         typename WaveFunctionTransfType,
         template<template<typename,typename,typename> class,
                  typename,
                  typename> class TargettingTemplate,
         typename MySparseMatrix>
void mainLoop2(GeometryType& geometry,
               const ParametersDmrgSolverType& dmrgSolverParams,
               InputNgType::Readable& io,
               const OperatorOptions& opOptions)
{
	if (dmrgSolverParams.options.find("ChebyshevSolver")!=PsimagLite::String::npos) {
		typedef TargettingTemplate<PsimagLite::ChebyshevSolver,
		                           MatrixVectorType,
		                           WaveFunctionTransfType> TargettingType;
		mainLoop3<GeometryType,TargettingType>
		(geometry,dmrgSolverParams,io,opOptions);
	} else {
		typedef TargettingTemplate<PsimagLite::LanczosSolver,
		                           MatrixVectorType,
		                           WaveFunctionTransfType> TargettingType;
		mainLoop3<GeometryType,TargettingType>
		(geometry,dmrgSolverParams,io,opOptions);
	}
}

template<typename GeometryType,
        template<typename> class ModelHelperTemplate,
         template<typename> class VectorWithOffsetTemplate,
         template<template<typename,typename,typename> class,
                  typename,
                  typename> class TargettingTemplate,
         typename MySparseMatrix>
void mainLoop(GeometryType& geometry,
              const ParametersDmrgSolverType& dmrgSolverParams,
              InputNgType::Readable& io,
              const OperatorOptions& opOptions)
{
	typedef Basis<MySparseMatrix> BasisType;
	typedef Operators<BasisType> OperatorsType;
	typedef BasisWithOperators<OperatorsType> BasisWithOperatorsType;
	typedef LeftRightSuper<BasisWithOperatorsType,BasisType> LeftRightSuperType;
	typedef ModelHelperTemplate<LeftRightSuperType> ModelHelperType;
	typedef ModelBase<ModelHelperType,
	                  ParametersDmrgSolverType,
	                  InputNgType::Readable,
	                  GeometryType> ModelBaseType;
	typedef typename BasisWithOperatorsType::SparseMatrixType SparseMatrixType;
	typedef typename SparseMatrixType::value_type ComplexOrRealType;
	typedef VectorWithOffsetTemplate<ComplexOrRealType> VectorWithOffsetType;
	typedef WaveFunctionTransfFactory<LeftRightSuperType,
	                                  VectorWithOffsetType> WaveFunctionTransfType;

	if (dmrgSolverParams.options.find("MatrixVectorStored")!=PsimagLite::String::npos) {
		mainLoop2<GeometryType,
		         MatrixVectorStored<ModelBaseType>,
		         WaveFunctionTransfType,
		         TargettingTemplate,
		         MySparseMatrix>(geometry,dmrgSolverParams,io,opOptions);
	} else if (dmrgSolverParams.options.find("MatrixVectorKron")!=PsimagLite::String::npos) {
		mainLoop2<GeometryType,
		MatrixVectorKron<ModelBaseType>,
		                 WaveFunctionTransfType,
		                 TargettingTemplate,
		                 MySparseMatrix>(geometry,dmrgSolverParams,io,opOptions);
	} else {
		mainLoop2<GeometryType,
		         MatrixVectorOnTheFly<ModelBaseType>,
		         WaveFunctionTransfType,
		         TargettingTemplate,
		         MySparseMatrix>(geometry,dmrgSolverParams,io,opOptions);
	}
}

template<typename GeometryType,
        template<typename> class ModelHelperTemplate,
         template<template<typename,typename,typename> class,
                  typename,
                  typename> class TargettingTemplate,
         typename MySparseMatrix>
void mainLoop1(GeometryType& geometry,
              const ParametersDmrgSolverType& dmrgSolverParams,
              InputNgType::Readable& io,
              const OperatorOptions& opOptions)
{
	if (dmrgSolverParams.options.find("vectorwithoffsets")!=PsimagLite::String::npos) {
		mainLoop<GeometryType,ModelHelperTemplate,VectorWithOffsets,TargettingTemplate,
	         MySparseMatrix>(geometry,dmrgSolverParams,io,opOptions);
	} else {
		mainLoop<GeometryType,ModelHelperTemplate,VectorWithOffset,TargettingTemplate,
	         MySparseMatrix>(geometry,dmrgSolverParams,io,opOptions);
	}
}

template<typename MySparseMatrix>
void mainLoop0(InputNgType::Readable& io,
               const ParametersDmrgSolverType& dmrgSolverParams,
               InputCheck& inputCheck,
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

	PsimagLite::String targetting=inputCheck.getTargeting(dmrgSolverParams.options);

	if (targetting!="GroundStateTargetting" && su2) {
		PsimagLite::String str("SU(2) supports only GroundStateTargetting (sorry!)\n");
		throw PsimagLite::RuntimeError(str);
	}

	if (su2) {
		mainLoop1<GeometryType,
		        ModelHelperSu2,
		        TargetingGroundState,
		        MySparseMatrix>(geometry,dmrgSolverParams,io,opOptions);

		return;
	}

#ifdef USE_COMPLEX
	if (targetting != "TimeStepTargetting" &&
	        targetting != "GroundStateTargetting") {
		PsimagLite::String str("USE_COMPLEX not allowed for ");
		str += targetting + "\n";
		throw PsimagLite::RuntimeError(str);
	}
#endif

	if (targetting=="TimeStepTargetting") {
		mainLoop1<GeometryType,ModelHelperLocal,TargetingTimeStep,
		         MySparseMatrixComplex>(geometry,dmrgSolverParams,io,opOptions);
		return;
	}

	if (targetting=="DynamicTargetting") {
		mainLoop1<GeometryType,ModelHelperLocal,TargetingDynamic,
		         MySparseMatrix>(geometry,dmrgSolverParams,io,opOptions);
		return;
	}

	if (targetting=="AdaptiveDynamicTargetting") {
		mainLoop1<GeometryType,ModelHelperLocal,TargetingAdaptiveDynamic,
		         MySparseMatrix>(geometry,dmrgSolverParams,io,opOptions);
		return;
	}

	if (targetting=="CorrectionVectorTargetting") {
		mainLoop1<GeometryType,ModelHelperLocal,TargetingCorrectionVector,
		         MySparseMatrix>(geometry,dmrgSolverParams,io,opOptions);
		return;
	}

	if (targetting=="CorrectionTargetting") {
		mainLoop1<GeometryType,ModelHelperLocal,TargetingCorrection,
		         MySparseMatrix>(geometry,dmrgSolverParams,io,opOptions);
		return;
	}

	if (targetting=="MettsTargetting") {
		mainLoop1<GeometryType,ModelHelperLocal,MettsTargetting,
		         MySparseMatrix>(geometry,dmrgSolverParams,io,opOptions);
		return;
	}

	if (targetting=="TargetingAncilla") {
		mainLoop1<GeometryType,ModelHelperLocal,TargetingTimeStep,
		         MySparseMatrix>(geometry,dmrgSolverParams,io,opOptions);
		return;
	}

	mainLoop1<GeometryType,ModelHelperLocal,TargetingGroundState,
	         MySparseMatrix>(geometry,dmrgSolverParams,io,opOptions);
}

int main(int argc,char *argv[])
{
	typedef PsimagLite::Concurrency ConcurrencyType;
	ConcurrencyType concurrency(&argc,&argv,1);
	InputCheck inputCheck;
	PsimagLite::String filename="";
	int opt = 0;
	OperatorOptions options;
	PsimagLite::String strUsage(argv[0]);
	if (utils::basename(strUsage) == "operator") options.enabled = true;
	strUsage += " -f filename";
	PsimagLite::String insitu("");
	int precision = 6;
	bool clearRun = false;
	/* PSIDOC DmrgDriver
	 * \begin{itemize}
	 * \item[-f] [Mandatory, String] Input to use.
	 * \item[-o] [Optional, String] What to measure in-situ
	 * \item[-l] [Optional, String] Redirect std::cout to a log file.
	 * If the string is ``.'' then choose name for log file automatically
	 * If the string is ``?'' then print name of log file that will
	 * be chosen if ``.'' is used.
	 *\item[-c] [Optional] Clear run for this input
	 * \end{itemize}
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
	./operator -l c -f input.inp -F -1\end{verbatim}
	See the function naturalOperator for each Model.

	\item[-d] [Optional, Integer] Degree of freedom (spin, orbital or
	combination of both) to use. This is model dependent. For example to
	obtain $c_\downarrow$ for the Hubbard model, say
	\begin{verbatim}./operator -l c -d 1 -f input.inp -F -1\end{verbatim}
	See the function naturalOperator for each Model. Defaults to 0.

	\item[-F] [Mandatory, 1 or -1] If this operator commutes on
	\emph{different} sites say 1 here. If it anticommutes on
	\emph{different} sites say -1.

	\item[-t] [Optional, Void] Transpose the operator. For example to
	obtain $c^\dagger_\uparrow$ for a Hubbard model, say
	\begin{verbatim}./operator -l c -t -f input.inp -F -1\end{verbatim}
	\end{itemize}
	 */
	while ((opt = getopt(argc, argv,"f:o:s:l:d:F:p:ct")) != -1) {
		switch (opt) {
		case 'f':
			filename = optarg;
			break;
		case 'o':
			if (insitu=="") {
				insitu = optarg;
			} else {
				insitu += ",";
				insitu += optarg;
			}
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
		case 'F':
			options.fermionicSign = atoi(optarg);
			break;
		case 'p':
			precision = atoi(optarg);
			std::cout.precision(precision);
			std::cerr.precision(precision);
			break;
		case 'c':
			clearRun = true;
			break;
		default:
			inputCheck.usageMain(strUsage);
			return 1;
		}
	}

	// sanity checks here
	if (filename=="") {
		inputCheck.usageMain(strUsage);
		return 1;
	}

	if (!options.enabled && options.label != "" && !clearRun) {
		bool queryOnly = (options.label == "?");
		if (options.label == "." || options.label == "?") {
			PsimagLite::String rootname = filename;
			size_t index =rootname.find(".", 0);
			if (index != PsimagLite::String::npos) {
				rootname.erase(index,filename.length());
			}

			options.label="runFor" + rootname + ".cout";
			if (queryOnly) {
				std::cout<<options.label<<"\n";
				return 0;
			}
		}

		GlobalCoutStream.open(options.label.c_str());
		if (!GlobalCoutStream || GlobalCoutStream.bad()
		        || !GlobalCoutStream.good()) {
			PsimagLite::String str(argv[0]);
			str += ": Could not redirect std::cout to " + options.label + "\n";
			throw PsimagLite::RuntimeError(str);
		}

		std::cerr<<argv[0]<<" ATTENTION: All standard output now sent to ";
		std::cerr<<options.label<<"\n";
		std::cerr.flush();
		GlobalCoutBuffer = std::cout.rdbuf(); //save old buf
		std::cout.rdbuf(GlobalCoutStream.rdbuf()); //redirect std::cout to file
		atexit(restoreCoutBuffer);
	}

	// print license
	if (ConcurrencyType::root()) {
		std::cout<<ProgramGlobals::license;
		Provenance provenance;
		std::cout<<provenance;
	}

	InputNgType::Writeable ioWriteable(filename,inputCheck);
	InputNgType::Readable io(ioWriteable);

	ParametersDmrgSolverType dmrgSolverParams(io);

	ArchiveFilesType af(dmrgSolverParams,filename,options.enabled,options.label);

	if (clearRun) {
		af.clear(filename);
		return 1;
	}

	if (insitu!="") dmrgSolverParams.insitu = insitu;

#ifndef USE_PTHREADS
	inputCheck.checkForThreads(dmrgSolverParams.nthreads);
#endif

	ConcurrencyType::npthreads = dmrgSolverParams.nthreads;

	registerSignals();

#ifdef USE_COMPLEX
	std::cerr<<argv[0]<<" EXPERIMENTAL option complex is in use\n";
	mainLoop0<MySparseMatrixC>(io,dmrgSolverParams,inputCheck,options);
#else
	mainLoop0<MySparseMatrixReal>(io,dmrgSolverParams,inputCheck,options);
#endif

	af.archiveIfNeeded();
}

