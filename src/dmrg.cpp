#include <string>
const std::string license=
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

typedef double MatrixElementType;

#include "CrsMatrix.h"
#include "LanczosSolver.h"
#include "ChebyshevSolver.h"
#include "BlockMatrix.h"
#include "DmrgSolver.h"
#include "IoSimple.h"
#include "Operator.h"
#ifndef USE_MPI
#include "ConcurrencySerial.h"
typedef PsimagLite::ConcurrencySerial<MatrixElementType> ConcurrencyType;
#else
#include "ConcurrencyMpi.h"
typedef PsimagLite::ConcurrencyMpi<MatrixElementType> ConcurrencyType;
#endif
#include "ModelFactory.h"
#include "OperatorsBase.h"
#include "Geometry.h"
#ifdef USE_PTHREADS
#include "Pthreads.h"
#define PTHREADS_NAME PsimagLite::Pthreads
#else
#include "NoPthreads.h"
#define PTHREADS_NAME PsimagLite::NoPthreads
#endif
#include "ModelHelperLocal.h"
#include "ModelHelperSu2.h"
#include "InternalProductOnTheFly.h"
#include "InternalProductStored.h"
#include "InternalProductCached.h"
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

typedef std::complex<MatrixElementType> ComplexType;
typedef  PsimagLite::CrsMatrix<ComplexType> MySparseMatrixComplex;
typedef  PsimagLite::CrsMatrix<MatrixElementType> MySparseMatrixReal;

using namespace Dmrg;

typedef PsimagLite::Geometry<MatrixElementType,ProgramGlobals> GeometryType;;
typedef PsimagLite::IoSimple IoType;
typedef IoType::In IoInputType;
typedef ParametersDmrgSolver<MatrixElementType> ParametersDmrgSolverType;

template<typename ModelFactoryType,
	 template<typename,typename> class InternalProductTemplate,
         typename TargettingType,
         typename MySparseMatrix>
void mainLoop3(GeometryType& geometry,
              ParametersDmrgSolverType& dmrgSolverParams,
              ConcurrencyType& concurrency,
              IoInputType& io)
{

	//! Setup the Model
	ModelFactoryType model(dmrgSolverParams,io,geometry);

	//! Read TimeEvolution if applicable:
	typedef typename TargettingType::TargettingParamsType TargettingParamsType;
	TargettingParamsType tsp(io,model);

	//! Setup the dmrg solver:
	typedef DmrgSolver<InternalProductTemplate,TargettingType> SolverType;
	SolverType dmrgSolver(dmrgSolverParams,model,concurrency,tsp);

	//! Calculate observables:
	dmrgSolver.main(geometry);
}

template<template<typename,typename> class ModelHelperTemplate,
	 template<typename,typename> class InternalProductTemplate,
         template<typename> class VectorWithOffsetTemplate,
         template<template<typename,typename,typename> class,
                  template<typename,typename> class,
                  template<typename,typename> class,
                  typename,typename,typename,
                  template<typename> class> class TargettingTemplate,
         typename MySparseMatrix>
void mainLoop2(GeometryType& geometry,
              ParametersDmrgSolverType& dmrgSolverParams,
              ConcurrencyType& concurrency,
              IoInputType& io)
{
	typedef Operator<MatrixElementType,MySparseMatrix> OperatorType;
	typedef Basis<MatrixElementType,MySparseMatrix> BasisType;
	typedef OperatorsBase<OperatorType,BasisType> OperatorsType;
	typedef BasisWithOperators<OperatorsType,ConcurrencyType> BasisWithOperatorsType;
	typedef LeftRightSuper<BasisWithOperatorsType,BasisType> LeftRightSuperType;
	typedef ModelHelperTemplate<LeftRightSuperType,ConcurrencyType> ModelHelperType;
	typedef ModelFactory<ModelHelperType,MySparseMatrix,GeometryType,
			     PTHREADS_NAME> ModelFactoryType;

	if (dmrgSolverParams.options.find("ChebyshevSolver")!=std::string::npos) {
		typedef TargettingTemplate<PsimagLite::ChebyshevSolver,
					   InternalProductTemplate,
					   WaveFunctionTransfFactory,
					   ModelFactoryType,
					   ConcurrencyType,
					   IoType,
					   VectorWithOffsetTemplate
					   > TargettingType;
		mainLoop3<ModelFactoryType, InternalProductTemplate,TargettingType,MySparseMatrix>
		(geometry,dmrgSolverParams,concurrency,io);
	} else {
		typedef TargettingTemplate<PsimagLite::LanczosSolver,
					   InternalProductTemplate,
					   WaveFunctionTransfFactory,
					   ModelFactoryType,
					   ConcurrencyType,
					   IoType,
					   VectorWithOffsetTemplate
					   > TargettingType;
		mainLoop3<ModelFactoryType,InternalProductTemplate,TargettingType,MySparseMatrix>
		(geometry,dmrgSolverParams,concurrency,io);
	}
}

template<template<typename,typename> class ModelHelperTemplate,
         template<typename> class VectorWithOffsetTemplate,
         template<template<typename,typename,typename> class,
                  template<typename,typename> class,
                  template<typename,typename> class,
                  typename,typename,typename,
                  template<typename> class> class TargettingTemplate,
         typename MySparseMatrix>
void mainLoop(GeometryType& geometry,
              ParametersDmrgSolverType& dmrgSolverParams,
              ConcurrencyType& concurrency,
              IoInputType& io)
{
	if (dmrgSolverParams.options.find("InternalProductStored")!=std::string::npos) {
		mainLoop2<ModelHelperTemplate,
		         InternalProductStored,
		         VectorWithOffsetTemplate,
		         TargettingTemplate,
		         MySparseMatrix>(geometry,dmrgSolverParams,concurrency,io);
	} else if (dmrgSolverParams.options.find("InternalProductCached")!=std::string::npos) {
		mainLoop2<ModelHelperTemplate,
		          InternalProductCached,
		          VectorWithOffsetTemplate,
		          TargettingTemplate,
		          MySparseMatrix>(geometry,dmrgSolverParams,concurrency,io);
	} else {
 		mainLoop2<ModelHelperTemplate,
		         InternalProductOnTheFly,
		         VectorWithOffsetTemplate,
		         TargettingTemplate,
		         MySparseMatrix>(geometry,dmrgSolverParams,concurrency,io);
	}
}

int main(int argc,char *argv[])
{
	//! setup distributed parallelization
	ConcurrencyType concurrency(argc,argv);

	// print license
	if (concurrency.root()) std::cerr<<license;

	//Setup the Geometry
	IoInputType io(argv[1]);
	GeometryType geometry(io);
	

	ParametersDmrgSolver<MatrixElementType> dmrgSolverParams(io);
	io.rewind();

	bool su2=false;
	if (dmrgSolverParams.options.find("useSu2Symmetry")!=std::string::npos) su2=true;
	std::string targetting="GroundStateTargetting";
	const char *targets[]={"TimeStepTargetting","DynamicTargetting","AdaptiveDynamicTargetting",
                     "CorrectionVectorTargetting","CorrectionTargetting","MettsTargetting"};
	size_t totalTargets = 6;
	for (size_t i = 0;i<totalTargets;++i)
		if (dmrgSolverParams.options.find(targets[i])!=std::string::npos) targetting=targets[i];

	if (targetting!="GroundStateTargetting" && su2) throw std::runtime_error("SU(2)"
 		" supports only GroundStateTargetting for now (sorry!)\n");
	if (su2) {
		if (dmrgSolverParams.targetQuantumNumbers[2]>0) { 
			mainLoop<ModelHelperSu2,VectorWithOffsets,GroundStateTargetting,
				MySparseMatrixReal>(geometry,dmrgSolverParams,concurrency,io);
		} else {
			mainLoop<ModelHelperSu2,VectorWithOffset,GroundStateTargetting,
				MySparseMatrixReal>(geometry,dmrgSolverParams,concurrency,io);
		}
		return 0;
	}
	if (targetting=="TimeStepTargetting") { 
		mainLoop<ModelHelperLocal,VectorWithOffsets,TimeStepTargetting,
			MySparseMatrixComplex>(geometry,dmrgSolverParams,concurrency,io);
			return 0;
	}
	if (targetting=="DynamicTargetting") {
		mainLoop<ModelHelperLocal,VectorWithOffsets,DynamicTargetting,
			MySparseMatrixReal>(geometry,dmrgSolverParams,concurrency,io);
			return 0;
	}
	if (targetting=="AdaptiveDynamicTargetting") {
		mainLoop<ModelHelperLocal,VectorWithOffsets,AdaptiveDynamicTargetting,
			MySparseMatrixReal>(geometry,dmrgSolverParams,concurrency,io);
			return 0;
	}
	if (targetting=="CorrectionVectorTargetting") {
		mainLoop<ModelHelperLocal,VectorWithOffsets,CorrectionVectorTargetting,
			MySparseMatrixReal>(geometry,dmrgSolverParams,concurrency,io);
			return 0;
	}
	if (targetting=="CorrectionTargetting") {
		mainLoop<ModelHelperLocal,VectorWithOffsets,CorrectionTargetting,
			MySparseMatrixReal>(geometry,dmrgSolverParams,concurrency,io);
			return 0;
	}
	if (targetting=="MettsTargetting") {
		mainLoop<ModelHelperLocal,VectorWithOffset,MettsTargetting,
			MySparseMatrixReal>(geometry,dmrgSolverParams,concurrency,io);
			return 0;
	}
	mainLoop<ModelHelperLocal,VectorWithOffset,GroundStateTargetting,
		MySparseMatrixReal>(geometry,dmrgSolverParams,concurrency,io);
}

