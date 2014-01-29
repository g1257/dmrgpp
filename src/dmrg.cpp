#include "String.h"

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

template<typename GeometryType,
         typename TargettingType>
void mainLoop3(GeometryType& geometry,
               const ParametersDmrgSolverType& dmrgSolverParams,
               InputNgType::Readable& io)
{
	typedef typename TargettingType::TargettingParamsType TargettingParamsType;
	typedef typename TargettingType::MatrixVectorType::ModelType ModelBaseType;

	//! Setup the Model
	ModelSelector<ModelBaseType> modelSelector(dmrgSolverParams.model);
	const ModelBaseType& model = modelSelector(dmrgSolverParams,io,geometry);

	//! Read TimeEvolution if applicable:
	TargettingParamsType tsp(io,model);

	//! Setup the dmrg solver:
	typedef DmrgSolver<TargettingType> SolverType;
	SolverType dmrgSolver(model,tsp);

	//! Calculate observables:
	dmrgSolver.main(geometry);
}

template<typename GeometryType,
         typename MatrixVectorType,
         typename WaveFunctionTransfType,
         template<template<typename,typename,typename> class,
                  typename,
                  typename,
                  typename> class TargettingTemplate,
         typename MySparseMatrix>
void mainLoop2(GeometryType& geometry,
               const ParametersDmrgSolverType& dmrgSolverParams,
               InputNgType::Readable& io)
{
	if (dmrgSolverParams.options.find("ChebyshevSolver")!=PsimagLite::String::npos) {
		typedef TargettingTemplate<PsimagLite::ChebyshevSolver,
					   MatrixVectorType,
					   WaveFunctionTransfType,
					   PsimagLite::IoSimple> TargettingType;
		mainLoop3<GeometryType,TargettingType>
		(geometry,dmrgSolverParams,io);
	} else {
		typedef TargettingTemplate<PsimagLite::LanczosSolver,
					   MatrixVectorType,
					   WaveFunctionTransfType,
					   PsimagLite::IoSimple> TargettingType;
		mainLoop3<GeometryType,TargettingType>
		(geometry,dmrgSolverParams,io);
	}
}

template<typename GeometryType,
        template<typename> class ModelHelperTemplate,
         template<typename> class VectorWithOffsetTemplate,
         template<template<typename,typename,typename> class,
                  typename,
                  typename,
                  typename> class TargettingTemplate,
         typename MySparseMatrix>
void mainLoop(GeometryType& geometry,
              const ParametersDmrgSolverType& dmrgSolverParams,
              InputNgType::Readable& io)
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
	typedef WaveFunctionTransfFactory<LeftRightSuperType,VectorWithOffsetType> WaveFunctionTransfType;

	if (dmrgSolverParams.options.find("MatrixVectorStored")!=PsimagLite::String::npos) {
		mainLoop2<GeometryType,
		         MatrixVectorStored<ModelBaseType>,
		         WaveFunctionTransfType,
		         TargettingTemplate,
		         MySparseMatrix>(geometry,dmrgSolverParams,io);
	} else if (dmrgSolverParams.options.find("MatrixVectorKron")!=PsimagLite::String::npos) {
		mainLoop2<GeometryType,
			 MatrixVectorKron<ModelBaseType>,
			 WaveFunctionTransfType,
			 TargettingTemplate,
			 MySparseMatrix>(geometry,dmrgSolverParams,io);
	} else {
 		mainLoop2<GeometryType,
		         MatrixVectorOnTheFly<ModelBaseType>,
		         WaveFunctionTransfType,
		         TargettingTemplate,
		         MySparseMatrix>(geometry,dmrgSolverParams,io);
	}
}


template<typename MySparseMatrix>
void mainLoop0(InputNgType::Readable& io,
               const ParametersDmrgSolverType& dmrgSolverParams,
               InputCheck& inputCheck)
{
	typedef typename MySparseMatrix::value_type ComplexOrRealType;
	typedef PsimagLite::Geometry<ComplexOrRealType,
	                             InputNgType::Readable,
	                             ProgramGlobals> GeometryType;

	GeometryType geometry(io);

	bool su2=false;
	if (dmrgSolverParams.options.find("useSu2Symmetry")!=PsimagLite::String::npos) su2=true;

	PsimagLite::String targetting=inputCheck.getTargeting(dmrgSolverParams.options);

	if (targetting!="GroundStateTargetting" && su2) throw PsimagLite::RuntimeError("SU(2)"
 		" supports only GroundStateTargetting for now (sorry!)\n");

	if (su2) {
		if (dmrgSolverParams.targetQuantumNumbers[2]>0) {
			mainLoop<GeometryType,ModelHelperSu2,VectorWithOffsets,TargetingGroundState,
				MySparseMatrix>(geometry,dmrgSolverParams,io);
		} else {
			mainLoop<GeometryType,ModelHelperSu2,VectorWithOffset,TargetingGroundState,
				MySparseMatrix>(geometry,dmrgSolverParams,io);
		}
		return;
	}

	if (targetting=="TimeStepTargetting") {
		mainLoop<GeometryType,ModelHelperLocal,VectorWithOffsets,TargetingTimeStep,
			MySparseMatrixComplex>(geometry,dmrgSolverParams,io);
			return;
	}
	if (targetting=="DynamicTargetting") {
		mainLoop<GeometryType,ModelHelperLocal,VectorWithOffsets,TargetingDynamic,
			MySparseMatrix>(geometry,dmrgSolverParams,io);
			return;
	}
	if (targetting=="AdaptiveDynamicTargetting") {
		mainLoop<GeometryType,ModelHelperLocal,VectorWithOffsets,TargetingAdaptiveDynamic,
			MySparseMatrix>(geometry,dmrgSolverParams,io);
			return;
	}
	if (targetting=="CorrectionVectorTargetting") {
		mainLoop<GeometryType,ModelHelperLocal,VectorWithOffsets,TargetingCorrectionVector,
			MySparseMatrix>(geometry,dmrgSolverParams,io);
			return;
	}
	if (targetting=="CorrectionTargetting") {
		mainLoop<GeometryType,ModelHelperLocal,VectorWithOffsets,TargetingCorrection,
			MySparseMatrix>(geometry,dmrgSolverParams,io);
			return;
	}
	if (targetting=="MettsTargetting") {
		mainLoop<GeometryType,ModelHelperLocal,VectorWithOffsets,MettsTargetting,
			MySparseMatrix>(geometry,dmrgSolverParams,io);
			return;
	}

	mainLoop<GeometryType,ModelHelperLocal,VectorWithOffset,TargetingGroundState,
		MySparseMatrix>(geometry,dmrgSolverParams,io);
}

int main(int argc,char *argv[])
{
	typedef PsimagLite::Concurrency ConcurrencyType;
	ConcurrencyType concurrency(&argc,&argv,1);
	InputCheck inputCheck;
	PsimagLite::String filename="";
	int opt = 0;
	PsimagLite::String strUsage(argv[0]);
	strUsage += " -f filename";
	PsimagLite::String insitu("");
	while ((opt = getopt(argc, argv,"f:o:")) != -1) {
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
		std::cout<<provenance;
	}

	InputNgType::Writeable ioWriteable(filename,inputCheck);
	InputNgType::Readable io(ioWriteable);

	ParametersDmrgSolverType dmrgSolverParams(io);

	if (insitu!="") dmrgSolverParams.insitu = insitu;

#ifndef USE_PTHREADS
	inputCheck.checkForThreads(dmrgSolverParams.nthreads);
#endif

	ConcurrencyType::npthreads = dmrgSolverParams.nthreads;

#ifdef USE_COMPLEX
		std::cerr<<argv[0]<<" EXPERIMENTAL option complex is in use\n";
		mainLoop0<MySparseMatrixC>(io,dmrgSolverParams,inputCheck);
#else
		mainLoop0<MySparseMatrixReal>(io,dmrgSolverParams,inputCheck);
#endif

}

