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
	while ((opt = getopt(argc, argv,"f:o:s:l:d:F:p:t")) != -1) {
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

	// print license
	if (ConcurrencyType::root()) {
		std::cerr<<ProgramGlobals::license;
		Provenance provenance;
		std::cerr<<provenance;
	}

	InputNgType::Writeable ioWriteable(filename,inputCheck);
	InputNgType::Readable io(ioWriteable);

	ParametersDmrgSolverType dmrgSolverParams(io);

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

}

