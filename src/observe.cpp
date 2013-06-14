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
#include "ModelFactory.h"
#include "OperatorsBase.h"
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
#include "MettsTargetting.h" // experimental
#include "BasisWithOperators.h"
#include "LeftRightSuper.h"
#include "InputNg.h"
#include "Provenance.h"

using namespace Dmrg;

typedef std::complex<RealType> ComplexType;

typedef  PsimagLite::CrsMatrix<ComplexType> MySparseMatrixComplex;
typedef  PsimagLite::CrsMatrix<RealType> MySparseMatrixReal;
typedef  PsimagLite::Geometry<RealType,ProgramGlobals> GeometryType;
typedef PsimagLite::IoSimple::In IoInputType;
typedef PsimagLite::InputNg<InputCheck> InputNgType;
typedef ParametersDmrgSolver<RealType,InputNgType::Readable> DmrgSolverParametersType;

template<typename ModelType>
SizeType dofsFromModelName(const ModelType& model)
{
	const PsimagLite::String& modelName = model.params().model;
	SizeType site = 0; // FIXME : account for Hilbert spaces changing with site
	SizeType dofs = SizeType(log(model.hilbertSize(site))/log(2.0));
	std::cerr<<"DOFS= "<<dofs<<" <------------------------------------\n";
	if (modelName.find("FeAsBasedSc")!=PsimagLite::String::npos) return dofs;
	if (modelName.find("FeAsBasedScExtended")!=PsimagLite::String::npos) return dofs;
	if (modelName.find("HubbardOneBand")!=PsimagLite::String::npos) return dofs;

	// max number here, site dependence taken into account elsewhere
	if (modelName.find("Immm")!=PsimagLite::String::npos) return 4;
	return 0;
}

template<typename VectorWithOffsetType,typename ModelType,typename SparseMatrixType,typename OperatorType,typename TargettingType>
bool observeOneFullSweep(
	IoInputType& io,
	const GeometryType& geometry,
	const ModelType& model,
	const PsimagLite::String& obsOptions,
	bool hasTimeEvolution)
{
	bool verbose = false;
	typedef typename SparseMatrixType::value_type FieldType;
	typedef Observer<FieldType,VectorWithOffsetType,ModelType,IoInputType> 
		ObserverType;
	typedef ObservableLibrary<ObserverType,TargettingType> ObservableLibraryType;
	SizeType n  = geometry.numberOfSites();
	//PsimagLite::String sSweeps = "sweeps=";
	//PsimagLite::String::SizeTypeype begin = obsOptions.find(sSweeps);
	//if (begin != PsimagLite::String::npos) {
	//	PsimagLite::String sTmp = obsOptions.substr(begin+sSweeps.length(),PsimagLite::String::npos);
		//std::cout<<"sTmp="<<sTmp<<"\n";
	//	n = atoi(sTmp.c_str());
	//}
	ObservableLibraryType observerLib(io,n,hasTimeEvolution,model,verbose);
	
	bool ot = false;
	if (obsOptions.find("ot")!=PsimagLite::String::npos || obsOptions.find("time")!=PsimagLite::String::npos) ot = true;
	if (hasTimeEvolution && ot) {
		observerLib.measureTime("superDensity");
		observerLib.measureTime("nupNdown");
		observerLib.measureTime("nup+ndown");
		if (obsOptions.find("sz")!=PsimagLite::String::npos) observerLib.measureTime("sz");
	}

	if (hasTimeEvolution) observerLib.setBrackets("time","time");

	const PsimagLite::String& modelName = model.params().model;
	SizeType rows = n; // could be n/2 if there's enough symmetry

	// Immm supports only onepoint:
	if (modelName=="Immm" && obsOptions!="onepoint") {
		PsimagLite::String str(__FILE__);
		str += " "  + ttos(__LINE__) + "\n";
		str += "Model Immm only supports onepoint\n";
		throw PsimagLite::RuntimeError(str.c_str());
	}

	SizeType numberOfDofs = dofsFromModelName(model);

	if (!hasTimeEvolution && obsOptions.find("onepoint")!=PsimagLite::String::npos) {
		observerLib.measureTheOnePoints(numberOfDofs);
	}

	if (modelName.find("Heisenberg")==PsimagLite::String::npos) {
		if (obsOptions.find("cc")!=PsimagLite::String::npos) {
			observerLib.measure("cc",rows,n);
		}

		if (obsOptions.find("nn")!=PsimagLite::String::npos) {
			observerLib.measure("nn",rows,n);
		}
	}
	if (obsOptions.find("szsz")!=PsimagLite::String::npos) {
		observerLib.measure("szsz",rows,n);
	}

	if (modelName.find("FeAsBasedSc")!=PsimagLite::String::npos ||
	    modelName.find("FeAsBasedScExtended")!=PsimagLite::String::npos ||
	    modelName.find("HubbardOneBand")!=PsimagLite::String::npos) {
		bool dd4 = (obsOptions.find("dd4")!=PsimagLite::String::npos);

		if (obsOptions.find("dd")!=PsimagLite::String::npos && !dd4) { // &&
			//geometry.label(0).find("ladder")!=PsimagLite::String::npos) {
			observerLib.measure("dd",rows,n);
		}

		// FOUR-POINT DELTA-DELTA^DAGGER:
		if (dd4 && geometry.label(0).find("ladder")!=PsimagLite::String::npos) {
			observerLib.measure("dd4",rows,n);
		} // if dd4
	}

	if (modelName.find("HubbardOneBand")!=PsimagLite::String::npos &&
	    obsOptions.find("multi")!=PsimagLite::String::npos) {
		observerLib.measure("multi",rows,n);
	}

	if (obsOptions.find("s+s-")!=PsimagLite::String::npos) {
		observerLib.measure("s+s-",rows,n);
	}

	if (obsOptions.find("ss")!=PsimagLite::String::npos) {
		observerLib.measure("ss",rows,n);
	}

	return observerLib.endOfData();
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
              const PsimagLite::String& targetting,
	      InputNgType::Readable& io,
	      const DmrgSolverParametersType& params,
              const PsimagLite::String& obsOptions)
{
	typedef Operator<RealType,MySparseMatrix> OperatorType;
	typedef Basis<RealType,MySparseMatrix> BasisType;
	typedef OperatorsBase<OperatorType,BasisType> OperatorsType;
	typedef typename OperatorType::SparseMatrixType SparseMatrixType;
	typedef BasisWithOperators<OperatorsType> BasisWithOperatorsType;
	typedef LeftRightSuper<BasisWithOperatorsType,BasisType> LeftRightSuperType;
	typedef ModelHelperTemplate<LeftRightSuperType> ModelHelperType;
	typedef ModelFactory<ModelHelperType,MySparseMatrix,GeometryType,PTHREADS_NAME,DmrgSolverParametersType> ModelType;
	typedef TargettingTemplate<PsimagLite::LanczosSolver,
	                           InternalProductOnTheFly,
	                           WaveFunctionTransfFactory,
	                           ModelType,
	                           IoInputType,
	                           VectorWithOffsetTemplate> TargettingType;

	typedef DmrgSolver<InternalProductOnTheFly,TargettingType> SolverType;
	
	typedef typename TargettingType::VectorWithOffsetType VectorWithOffsetType;
	
	ModelType model(params,io,geometry);

	 //! Read TimeEvolution if applicable:
	typedef typename TargettingType::TargettingParamsType TargettingParamsType;
	TargettingParamsType tsp(io,model);
	
	bool moreData = true;
	const PsimagLite::String& datafile = params.filename;
	IoInputType dataIo(datafile);
	bool hasTimeEvolution = (targetting == "TimeStepTargetting" || targetting=="MettsTargetting") ? true : false;
	while (moreData) {
		try {
			moreData = !observeOneFullSweep<VectorWithOffsetType,ModelType,
			            SparseMatrixType,OperatorType,TargettingType>
			(dataIo,geometry,model,obsOptions,hasTimeEvolution);
		} catch (std::exception& e) {
			std::cerr<<"CAUGHT: "<<e.what();
			std::cerr<<"There's no more data\n";
			break;
		}

		//if (!hasTimeEvolution) break;
	}
}

void usage(const char* name)
{
	std::cerr<<"USAGE is "<<name<<" -f filename [-o options]\n";
}

int main(int argc,char *argv[])
{
	using namespace Dmrg;

	PsimagLite::String filename="";
	PsimagLite::String options = "";
	int opt = 0;
	while ((opt = getopt(argc, argv,"f:o:")) != -1) {
		switch (opt) {
		case 'f':
			filename = optarg;
			break;
		case 'o':
			options = optarg;
			break;
		default:
			usage(argv[0]);
			return 1;
		}
	}

	//sanity checks here
	if (filename=="") {
		usage(argv[0]);
		return 1;
	}

	typedef PsimagLite::Concurrency ConcurrencyType;
	ConcurrencyType concurrency(argc,argv);
	
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
	if (targetting=="MettsTargetting") { // experimental, do not use
		mainLoop<ModelHelperLocal,VectorWithOffsets,MettsTargetting,MySparseMatrixReal>
		(geometry,targetting,io,dmrgSolverParams,options);
		return 0;
	}
	mainLoop<ModelHelperLocal,VectorWithOffset,GroundStateTargetting,MySparseMatrixReal>
	(geometry,targetting,io,dmrgSolverParams,options);
} // main

