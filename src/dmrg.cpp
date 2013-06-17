#include "String.h"
const PsimagLite::String license=
"Copyright (c) 2009-2012 , UT-Battelle, LLC\n"
"All rights reserved\n"
"\n"
"[DMRG++, Version 2.0.0]\n"
"\n"
"*********************************************************\n"
"THE SOFTWARE IS SUPPLIED BY THE COPYRIGHT HOLDERS AND\n"
"CONTRIBUTORS \"AS IS\" AND ANY EXPRESS OR IMPLIED\n"
"WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED\n"
"WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A\n"
"PARTICULAR PURPOSE ARE DISCLAIMED.\n"
"\n"
"Please see full open source license included in file LICENSE.\n"
"*********************************************************\n"
"\n";

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
#include "Operator.h"
#include "ModelFactory.h"
#include "OperatorsBase.h"
#include "Concurrency.h"
#include "Geometry/Geometry.h"
#include "ModelHelperLocal.h"
#include "ModelHelperSu2.h"
#include "InternalProductOnTheFly.h"
#include "InternalProductStored.h"
#include "InternalProductKron.h"
//#include "InternalProductCached.h"
#include "GroundStateTargetting.h"
#include "TimeStepTargetting.h"
#include "DynamicTargetting.h"
#include "AdaptiveDynamicTargetting.h"
#include "CorrectionTargetting.h"
#include "CorrectionVectorTargetting.h"
#include "MettsTargetting.h" // experimental
#include "VectorWithOffset.h"
#include "VectorWithOffsets.h"
#include "BasisWithOperators.h"
#include "LeftRightSuper.h"
#include "Provenance.h"

typedef std::complex<MatrixElementType> ComplexType;
typedef  PsimagLite::CrsMatrix<ComplexType> MySparseMatrixComplex;
typedef  PsimagLite::CrsMatrix<MatrixElementType> MySparseMatrixReal;

using namespace Dmrg;

typedef PsimagLite::Geometry<MatrixElementType,ProgramGlobals> GeometryType;
typedef PsimagLite::InputNg<InputCheck> InputNgType;
typedef ParametersDmrgSolver<MatrixElementType,InputNgType::Readable> ParametersDmrgSolverType;

template<typename ModelFactoryType,
	 template<typename,typename> class InternalProductTemplate,
         typename TargettingType>
void mainLoop3(GeometryType& geometry,
               ParametersDmrgSolverType& dmrgSolverParams,
               InputNgType::Readable& io)
{
	//! Setup the Model
	ModelFactoryType model(dmrgSolverParams,io,geometry);

	//! Read TimeEvolution if applicable:
	typedef typename TargettingType::TargettingParamsType TargettingParamsType;
	TargettingParamsType tsp(io,model);

	//! Setup the dmrg solver:
	typedef DmrgSolver<InternalProductTemplate,TargettingType> SolverType;
	SolverType dmrgSolver(dmrgSolverParams,model,tsp);

	//! Calculate observables:
	dmrgSolver.main(geometry);
}

template<template<typename> class ModelHelperTemplate,
	 template<typename,typename> class InternalProductTemplate,
         template<typename> class VectorWithOffsetTemplate,
         template<template<typename,typename,typename> class,
                  template<typename,typename> class,
                  template<typename,typename> class,
                  typename,typename,
                  template<typename> class> class TargettingTemplate,
         typename MySparseMatrix>
void mainLoop2(GeometryType& geometry,
               ParametersDmrgSolverType& dmrgSolverParams,
               InputNgType::Readable& io)
{
	typedef Operator<MatrixElementType,MySparseMatrix> OperatorType;
	typedef Basis<MatrixElementType,MySparseMatrix> BasisType;
	typedef OperatorsBase<OperatorType,BasisType> OperatorsType;
	typedef BasisWithOperators<OperatorsType> BasisWithOperatorsType;
	typedef LeftRightSuper<BasisWithOperatorsType,BasisType> LeftRightSuperType;
	typedef ModelHelperTemplate<LeftRightSuperType> ModelHelperType;
	typedef ModelFactory<ModelHelperType,MySparseMatrix,GeometryType,
	        ParametersDmrgSolverType> ModelFactoryType;

	if (dmrgSolverParams.options.find("ChebyshevSolver")!=PsimagLite::String::npos) {
		typedef TargettingTemplate<PsimagLite::ChebyshevSolver,
					   InternalProductTemplate,
					   WaveFunctionTransfFactory,
					   ModelFactoryType,
					   PsimagLite::IoSimple,
					   VectorWithOffsetTemplate
					   > TargettingType;
		mainLoop3<ModelFactoryType, InternalProductTemplate,TargettingType>
		(geometry,dmrgSolverParams,io);
	} else {
		typedef TargettingTemplate<PsimagLite::LanczosSolver,
					   InternalProductTemplate,
					   WaveFunctionTransfFactory,
					   ModelFactoryType,
					   PsimagLite::IoSimple,
					   VectorWithOffsetTemplate
					   > TargettingType;
		mainLoop3<ModelFactoryType,InternalProductTemplate,TargettingType>
		(geometry,dmrgSolverParams,io);
	}
}

template<template<typename> class ModelHelperTemplate,
         template<typename> class VectorWithOffsetTemplate,
         template<template<typename,typename,typename> class,
                  template<typename,typename> class,
                  template<typename,typename> class,
                  typename,typename,
                  template<typename> class> class TargettingTemplate,
         typename MySparseMatrix>
void mainLoop(GeometryType& geometry,
              ParametersDmrgSolverType& dmrgSolverParams,
	      InputNgType::Readable& io)
{
	if (dmrgSolverParams.options.find("InternalProductStored")!=PsimagLite::String::npos) {
		mainLoop2<ModelHelperTemplate,
		         InternalProductStored,
		         VectorWithOffsetTemplate,
		         TargettingTemplate,
		         MySparseMatrix>(geometry,dmrgSolverParams,io);
	} else if (dmrgSolverParams.options.find("InternalProductKron")!=PsimagLite::String::npos) {
		mainLoop2<ModelHelperTemplate,
			 InternalProductKron,
			 VectorWithOffsetTemplate,
			 TargettingTemplate,
			 MySparseMatrix>(geometry,dmrgSolverParams,io);
	} else {
 		mainLoop2<ModelHelperTemplate,
		         InternalProductOnTheFly,
		         VectorWithOffsetTemplate,
		         TargettingTemplate,
		         MySparseMatrix>(geometry,dmrgSolverParams,io);
	}
}

int main(int argc,char *argv[])
{
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

	typedef PsimagLite::Concurrency ConcurrencyType;
	ConcurrencyType concurrency(argc,argv,1);

	// print license
	if (ConcurrencyType::root()) {
		std::cerr<<license;
		Provenance provenance;
		std::cout<<provenance;
	}

	InputNgType::Writeable ioWriteable(filename,inputCheck);
	InputNgType::Readable io(ioWriteable);

	GeometryType geometry(io);

	ParametersDmrgSolver<MatrixElementType,InputNgType::Readable> dmrgSolverParams(io);

	if (insitu!="") dmrgSolverParams.insitu = insitu;

#ifndef USE_PTHREADS
	inputCheck.checkForThreads(dmrgSolverParams.nthreads);
#endif

	ConcurrencyType::npthreads = dmrgSolverParams.nthreads;

	bool su2=false;
	if (dmrgSolverParams.options.find("useSu2Symmetry")!=PsimagLite::String::npos) su2=true;
	PsimagLite::String targetting="GroundStateTargetting";
	const char *targets[]={"TimeStepTargetting","DynamicTargetting","AdaptiveDynamicTargetting",
                     "CorrectionVectorTargetting","CorrectionTargetting","MettsTargetting"};
	SizeType totalTargets = 6;
	for (SizeType i = 0;i<totalTargets;++i)
		if (dmrgSolverParams.options.find(targets[i])!=PsimagLite::String::npos) targetting=targets[i];

	if (targetting!="GroundStateTargetting" && su2) throw PsimagLite::RuntimeError("SU(2)"
 		" supports only GroundStateTargetting for now (sorry!)\n");
	if (su2) {
		if (dmrgSolverParams.targetQuantumNumbers[2]>0) { 
			mainLoop<ModelHelperSu2,VectorWithOffsets,GroundStateTargetting,
				MySparseMatrixReal>(geometry,dmrgSolverParams,io);
		} else {
			mainLoop<ModelHelperSu2,VectorWithOffset,GroundStateTargetting,
				MySparseMatrixReal>(geometry,dmrgSolverParams,io);
		}
		return 0;
	}
	if (targetting=="TimeStepTargetting") { 
		mainLoop<ModelHelperLocal,VectorWithOffsets,TimeStepTargetting,
			MySparseMatrixComplex>(geometry,dmrgSolverParams,io);
			return 0;
	}
	if (targetting=="DynamicTargetting") {
		mainLoop<ModelHelperLocal,VectorWithOffsets,DynamicTargetting,
			MySparseMatrixReal>(geometry,dmrgSolverParams,io);
			return 0;
	}
	if (targetting=="AdaptiveDynamicTargetting") {
		mainLoop<ModelHelperLocal,VectorWithOffsets,AdaptiveDynamicTargetting,
			MySparseMatrixReal>(geometry,dmrgSolverParams,io);
			return 0;
	}
	if (targetting=="CorrectionVectorTargetting") {
		mainLoop<ModelHelperLocal,VectorWithOffsets,CorrectionVectorTargetting,
			MySparseMatrixReal>(geometry,dmrgSolverParams,io);
			return 0;
	}
	if (targetting=="CorrectionTargetting") {
		mainLoop<ModelHelperLocal,VectorWithOffsets,CorrectionTargetting,
			MySparseMatrixReal>(geometry,dmrgSolverParams,io);
			return 0;
	}
	if (targetting=="MettsTargetting") {
		mainLoop<ModelHelperLocal,VectorWithOffsets,MettsTargetting,
			MySparseMatrixReal>(geometry,dmrgSolverParams,io);
			return 0;
	}
	mainLoop<ModelHelperLocal,VectorWithOffset,GroundStateTargetting,
		MySparseMatrixReal>(geometry,dmrgSolverParams,io);
}

