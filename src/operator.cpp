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
typedef double RealType;
#else
typedef float RealType;
#endif

#include <unistd.h>
#include "Observer.h"
#include "ObservableLibrary.h"
#include "IoSimple.h"
#include "Operators.h"
#include "Concurrency.h"
#include "Geometry/Geometry.h" 
#include "CrsMatrix.h"
#include "ModelHelperLocal.h"
#include "ModelHelperSu2.h"
#include "InternalProductOnTheFly.h"
#include "VectorWithOffset.h"
#include "VectorWithOffsets.h"
#include "GroundStateTargetting.h"
#include "DmrgSolver.h" // only used for types
#include "TimeStepTargetting.h"
#include "DynamicTargetting.h"
#include "AdaptiveDynamicTargetting.h"
#include "CorrectionTargetting.h"
#include "MettsTargetting.h"
#include "BasisWithOperators.h"
#include "LeftRightSuper.h"
#include "InputNg.h"
#include "Provenance.h"
#include "InputCheck.h"
#include "ModelSelector.h"

using namespace Dmrg;

typedef std::complex<RealType> ComplexType;

typedef  PsimagLite::CrsMatrix<ComplexType> MySparseMatrixComplex;
typedef  PsimagLite::CrsMatrix<RealType> MySparseMatrixReal;
typedef  PsimagLite::Geometry<RealType,ProgramGlobals> GeometryType;
typedef PsimagLite::IoSimple::In IoInputType;
typedef PsimagLite::InputNg<InputCheck> InputNgType;
typedef ParametersDmrgSolver<RealType,InputNgType::Readable> DmrgSolverParametersType;

struct OperatorOptions {
	
	OperatorOptions()
	: site(0),dof(0),label("")
	{}
	
	SizeType site;
	SizeType dof;
	PsimagLite::String label;
};

template<template<typename> class ModelHelperTemplate,
         template<typename> class VectorWithOffsetTemplate,
         template<template<typename,typename,typename> class,
                  template<typename,typename> class,
                  template<typename,typename> class,
                  typename,typename,
                  template<typename> class> class TargettingTemplate,
         typename MySparseMatrix>
void mainLoop(GeometryType& geometry,
              const PsimagLite::String& targetting,
              InputNgType::Readable& io,
              const DmrgSolverParametersType& params,
              const OperatorOptions& obsOptions)
{
	typedef Basis<MySparseMatrix> BasisType;
	typedef Operators<BasisType> OperatorsType;
	typedef BasisWithOperators<OperatorsType> BasisWithOperatorsType;
	typedef LeftRightSuper<BasisWithOperatorsType,BasisType> LeftRightSuperType;
	typedef ModelHelperTemplate<LeftRightSuperType> ModelHelperType;
	typedef ModelBase<ModelHelperType,
	                  DmrgSolverParametersType,
	                  InputNgType::Readable,
	                  GeometryType> ModelBaseType;
	typedef TargettingTemplate<PsimagLite::LanczosSolver,
	                           InternalProductOnTheFly,
	                           WaveFunctionTransfFactory,
	                           ModelBaseType,
	                           IoInputType,
	                           VectorWithOffsetTemplate> TargettingType;

	typedef DmrgSolver<InternalProductOnTheFly,TargettingType> SolverType;
	
	typedef typename TargettingType::VectorWithOffsetType VectorWithOffsetType;
	typedef typename ModelHelperType::SparseElementType SparseElementType;
	typedef PsimagLite::Matrix<SparseElementType> MatrixType;
	
	ModelSelector<ModelBaseType> modelSelector(params.model);
	const ModelBaseType& model = modelSelector(params,io,geometry);
	

	MatrixType opC2 = model.naturalOperator(obsOptions.label,obsOptions.site,obsOptions.dof);
	std::cout<<"label="<<obsOptions.label<<" site="<<obsOptions.site<<" dof="<<obsOptions.dof<<"\n";
	std::cout<<opC2;
}

void usage(const char* name)
{
	std::cerr<<"USAGE is "<<name<<" -f filename -s site -l label -d dof\n";
}

int main(int argc,char *argv[])
{
	typedef PsimagLite::Concurrency ConcurrencyType;

	using namespace Dmrg;

	PsimagLite::String filename="";
	OperatorOptions options;
	int opt = 0;
	while ((opt = getopt(argc, argv,"f:s:l:d:")) != -1) {
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
		default:
			usage(argv[0]);
			return 1;
		}
	}

	//sanity checks here
	if (filename=="" || options.label=="") {
		usage(argv[0]);
		return 1;
	}

	ConcurrencyType concurrency(&argc,&argv,1);
	
	// print license
	if (ConcurrencyType::root()) {
		std::cerr<<license;
		Provenance provenance;
		std::cout<<provenance;
	}

	//Setup the Geometry
	InputCheck inputCheck;
	InputNgType::Writeable ioWriteable(filename,inputCheck);
	InputNgType::Readable io(ioWriteable);
	GeometryType geometry(io);

	//! Read the parameters for this run
	//ParametersModelType mp(io);
	DmrgSolverParametersType dmrgSolverParams(io);


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
			mainLoop<ModelHelperSu2,VectorWithOffsets,GroundStateTargetting,MySparseMatrixReal>
			(geometry,targetting,io,dmrgSolverParams,options);
		} else {
			mainLoop<ModelHelperSu2,VectorWithOffset,GroundStateTargetting,MySparseMatrixReal>
			(geometry,targetting,io,dmrgSolverParams,options);
		}
		return 0;
	}
	if (targetting=="TimeStepTargetting") { 
		mainLoop<ModelHelperLocal,VectorWithOffsets,TimeStepTargetting,MySparseMatrixComplex>
		(geometry,targetting,io,dmrgSolverParams,options);
		return 0;
	}
	if (targetting=="DynamicTargetting") {
		mainLoop<ModelHelperLocal,VectorWithOffsets,DynamicTargetting,MySparseMatrixReal>
		(geometry,targetting,io,dmrgSolverParams,options);
		return 0;
	}
	if (targetting=="AdaptiveDynamicTargetting") {
		mainLoop<ModelHelperLocal,VectorWithOffsets,AdaptiveDynamicTargetting,MySparseMatrixReal>
		(geometry,targetting,io,dmrgSolverParams,options);
		return 0;
	}
	if (targetting=="CorrectionTargetting") {
		mainLoop<ModelHelperLocal,VectorWithOffsets,CorrectionTargetting,MySparseMatrixReal>
		(geometry,targetting,io,dmrgSolverParams,options);
		return 0;
	}
	if (targetting=="MettsTargetting") {
		mainLoop<ModelHelperLocal,VectorWithOffset,MettsTargetting,MySparseMatrixReal>
		(geometry,targetting,io,dmrgSolverParams,options);
		return 0;
	}
	mainLoop<ModelHelperLocal,VectorWithOffset,GroundStateTargetting,MySparseMatrixReal>
	(geometry,targetting,io,dmrgSolverParams,options);
} // main

