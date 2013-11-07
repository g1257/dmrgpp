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
#include "Operators.h"
#include "Concurrency.h"
#include "Geometry/Geometry.h"
#include "ModelHelperLocal.h"
#include "ModelHelperSu2.h"
#include "InternalProductOnTheFly.h"
#include "InternalProductStored.h"
#include "InternalProductKron.h"
#include "GroundStateTargetting.h"
#include "TimeStepTargetting.h"
#include "DynamicTargetting.h"
#include "AdaptiveDynamicTargetting.h"
#include "CorrectionTargetting.h"
#include "CorrectionVectorTargetting.h"
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
         typename ModelBaseType,
         template<typename,typename> class InternalProductTemplate,
         typename TargettingType>
void mainLoop3(GeometryType& geometry,
               ParametersDmrgSolverType& dmrgSolverParams,
               InputNgType::Readable& io)
{
	//! Setup the Model
	ModelSelector<ModelBaseType> modelSelector(dmrgSolverParams.model);
	const ModelBaseType& model = modelSelector(dmrgSolverParams,io,geometry);

	//! Read TimeEvolution if applicable:
	typedef typename TargettingType::TargettingParamsType TargettingParamsType;
	TargettingParamsType tsp(io,model);

	//! Setup the dmrg solver:
	typedef DmrgSolver<InternalProductTemplate,TargettingType> SolverType;
	SolverType dmrgSolver(dmrgSolverParams,model,tsp);

	//! Calculate observables:
	dmrgSolver.main(geometry);
}

template<typename GeometryType,
         template<typename> class ModelHelperTemplate,
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
	typedef Basis<MySparseMatrix> BasisType;
	typedef Operators<BasisType> OperatorsType;
	typedef BasisWithOperators<OperatorsType> BasisWithOperatorsType;
	typedef LeftRightSuper<BasisWithOperatorsType,BasisType> LeftRightSuperType;
	typedef ModelHelperTemplate<LeftRightSuperType> ModelHelperType;
	typedef ModelBase<ModelHelperType,
	                  ParametersDmrgSolverType,
	                  InputNgType::Readable,
	                  GeometryType> ModelBaseType;

	if (dmrgSolverParams.options.find("ChebyshevSolver")!=PsimagLite::String::npos) {
		typedef TargettingTemplate<PsimagLite::ChebyshevSolver,
					   InternalProductTemplate,
					   WaveFunctionTransfFactory,
					   ModelBaseType,
					   PsimagLite::IoSimple,
					   VectorWithOffsetTemplate
					   > TargettingType;
		mainLoop3<GeometryType,ModelBaseType,InternalProductTemplate,TargettingType>
		(geometry,dmrgSolverParams,io);
	} else {
		typedef TargettingTemplate<PsimagLite::LanczosSolver,
					   InternalProductTemplate,
					   WaveFunctionTransfFactory,
					   ModelBaseType,
					   PsimagLite::IoSimple,
					   VectorWithOffsetTemplate
					   > TargettingType;
		mainLoop3<GeometryType,ModelBaseType,InternalProductTemplate,TargettingType>
		(geometry,dmrgSolverParams,io);
	}
}

template<typename GeometryType,
        template<typename> class ModelHelperTemplate,
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
		mainLoop2<GeometryType,
		         ModelHelperTemplate,
		         InternalProductStored,
		         VectorWithOffsetTemplate,
		         TargettingTemplate,
		         MySparseMatrix>(geometry,dmrgSolverParams,io);
	} else if (dmrgSolverParams.options.find("InternalProductKron")!=PsimagLite::String::npos) {
		mainLoop2<GeometryType,
		     ModelHelperTemplate,
			 InternalProductKron,
			 VectorWithOffsetTemplate,
			 TargettingTemplate,
			 MySparseMatrix>(geometry,dmrgSolverParams,io);
	} else {
 		mainLoop2<GeometryType,
		         ModelHelperTemplate,
		         InternalProductOnTheFly,
		         VectorWithOffsetTemplate,
		         TargettingTemplate,
		         MySparseMatrix>(geometry,dmrgSolverParams,io);
	}
}


template<typename MySparseMatrix>
void mainLoop0(InputNgType::Readable& io,
               ParametersDmrgSolverType& dmrgSolverParams,
               InputCheck& inputCheck)
{
	typedef typename MySparseMatrix::value_type ComplexOrRealType;
	typedef PsimagLite::Geometry<ComplexOrRealType,ProgramGlobals> GeometryType;

	GeometryType geometry(io);

	bool su2=false;
	if (dmrgSolverParams.options.find("useSu2Symmetry")!=PsimagLite::String::npos) su2=true;

	PsimagLite::String targetting=inputCheck.getTargeting(dmrgSolverParams.options);

	if (targetting!="GroundStateTargetting" && su2) throw PsimagLite::RuntimeError("SU(2)"
 		" supports only GroundStateTargetting for now (sorry!)\n");

	if (su2) {
		if (dmrgSolverParams.targetQuantumNumbers[2]>0) {
			mainLoop<GeometryType,ModelHelperSu2,VectorWithOffsets,GroundStateTargetting,
				MySparseMatrix>(geometry,dmrgSolverParams,io);
		} else {
			mainLoop<GeometryType,ModelHelperSu2,VectorWithOffset,GroundStateTargetting,
				MySparseMatrix>(geometry,dmrgSolverParams,io);
		}
		return;
	}

	if (targetting=="TimeStepTargetting") {
		mainLoop<GeometryType,ModelHelperLocal,VectorWithOffsets,TimeStepTargetting,
			MySparseMatrixComplex>(geometry,dmrgSolverParams,io);
			return;
	}
	if (targetting=="DynamicTargetting") {
		mainLoop<GeometryType,ModelHelperLocal,VectorWithOffsets,DynamicTargetting,
			MySparseMatrix>(geometry,dmrgSolverParams,io);
			return;
	}
	if (targetting=="AdaptiveDynamicTargetting") {
		mainLoop<GeometryType,ModelHelperLocal,VectorWithOffsets,AdaptiveDynamicTargetting,
			MySparseMatrix>(geometry,dmrgSolverParams,io);
			return;
	}
	if (targetting=="CorrectionVectorTargetting") {
		mainLoop<GeometryType,ModelHelperLocal,VectorWithOffsets,CorrectionVectorTargetting,
			MySparseMatrix>(geometry,dmrgSolverParams,io);
			return;
	}
	if (targetting=="CorrectionTargetting") {
		mainLoop<GeometryType,ModelHelperLocal,VectorWithOffsets,CorrectionTargetting,
			MySparseMatrix>(geometry,dmrgSolverParams,io);
			return;
	}
	if (targetting=="MettsTargetting") {
		mainLoop<GeometryType,ModelHelperLocal,VectorWithOffsets,MettsTargetting,
			MySparseMatrix>(geometry,dmrgSolverParams,io);
			return;
	}

	mainLoop<GeometryType,ModelHelperLocal,VectorWithOffset,GroundStateTargetting,
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
		std::cerr<<license;
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

