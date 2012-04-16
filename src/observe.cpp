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

typedef double RealType;

#include "Observer.h"
#include "ObservableLibrary.h"
#include "IoSimple.h"
#include "ModelFactory.h"
#include "OperatorsBase.h"
#ifndef USE_MPI
#include "ConcurrencySerial.h"
typedef PsimagLite::ConcurrencySerial<RealType> ConcurrencyType;
#else
#include "ConcurrencyMpi.h"
typedef PsimagLite::ConcurrencyMpi<RealType> ConcurrencyType;
#endif
#include "Geometry.h" 
#include "CrsMatrix.h"
#ifdef USE_PTHREADS
#include "Pthreads.h"
#define PTHREADS_NAME Pthreads
#else
#include "NoPthreads.h"
#define PTHREADS_NAME NoPthreads
#endif 
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
#include "InputValidator.h"
#include "LabelWithKnownSize.h"

using namespace Dmrg;

typedef std::complex<RealType> ComplexType;

typedef  PsimagLite::CrsMatrix<ComplexType> MySparseMatrixComplex;
typedef  PsimagLite::CrsMatrix<RealType> MySparseMatrixReal;
typedef  PsimagLite::Geometry<RealType,ProgramGlobals> GeometryType;
typedef PsimagLite::IoSimple::In IoInputType;
typedef ParametersDmrgSolver<RealType> DmrgSolverParametersType;
typedef PsimagLite::LabelWithKnownSize LabelWithKnownSizeType;

size_t dofsFromModelName(const std::string& modelName)
{
	if (modelName.find("FeAsBasedSc")!=std::string::npos) return 4;
	if (modelName.find("FeAsBasedScExtended")!=std::string::npos) return 4;
	if (modelName.find("HubbardOneBand")!=std::string::npos) return 2;
	return 0;
}

template<typename VectorWithOffsetType,typename ModelType,typename SparseMatrixType,typename OperatorType,typename TargettingType>
bool observeOneFullSweep(
	IoInputType& io,
	const GeometryType& geometry,
	const ModelType& model,
	const std::string& obsOptions,
	bool hasTimeEvolution,
	ConcurrencyType& concurrency)
{
	bool verbose = false;
	typedef typename SparseMatrixType::value_type FieldType;
	typedef Observer<FieldType,VectorWithOffsetType,ModelType,IoInputType> 
		ObserverType;
	typedef ObservableLibrary<ObserverType,TargettingType> ObservableLibraryType;
	size_t n  = geometry.numberOfSites();
	//std::string sSweeps = "sweeps=";
	//std::string::size_type begin = obsOptions.find(sSweeps);
	//if (begin != std::string::npos) {
	//	std::string sTmp = obsOptions.substr(begin+sSweeps.length(),std::string::npos);
		//std::cout<<"sTmp="<<sTmp<<"\n";
	//	n = atoi(sTmp.c_str());
	//}
	ObservableLibraryType observerLib(io,n,hasTimeEvolution,model,concurrency,verbose);
	
	bool ot = false;
	if (obsOptions.find("ot")!=std::string::npos || obsOptions.find("time")!=std::string::npos) ot = true;
	if (hasTimeEvolution && ot) {
		observerLib.measureTime("superDensity");
		observerLib.measureTime("nupNdown");
		observerLib.measureTime("nup+ndown");
		if (obsOptions.find("sz")!=std::string::npos) observerLib.measureTime("sz");
	}

	if (hasTimeEvolution) observerLib.setBrackets("time","time");

	const std::string& modelName = model.name();
	size_t rows = n; // could be n/2 if there's enough symmetry

	size_t numberOfDofs = dofsFromModelName(modelName);
	if (!hasTimeEvolution && obsOptions.find("onepoint")!=std::string::npos) {
		observerLib.measureTheOnePoints(numberOfDofs);
	}

	if (modelName.find("Heisenberg")==std::string::npos) {
		if (obsOptions.find("cc")!=std::string::npos) {
			observerLib.measure("cc",rows,n);
		}

		if (obsOptions.find("nn")!=std::string::npos) {
			observerLib.measure("nn",rows,n);
		}
	}
	if (obsOptions.find("szsz")!=std::string::npos) {
		observerLib.measure("szsz",rows,n);
	}

	if (modelName.find("FeAsBasedSc")!=std::string::npos ||
	    modelName.find("FeAsBasedScExtended")!=std::string::npos) {
		if (obsOptions.find("dd")!=std::string::npos &&
		    geometry.label(0).find("ladder")!=std::string::npos) {
			observerLib.measure("dd",rows,n);
		}

		// FOUR-POINT DELTA-DELTA^DAGGER:
		if (obsOptions.find("dd4")!=std::string::npos &&
		    geometry.label(0).find("ladder")!=std::string::npos) {
			observerLib.measure("dd4",rows,n);
		} // if dd4
	}

	//if (s.find("heisenberg")!=std::string::npos) {
	if (obsOptions.find("s+s-")!=std::string::npos) {
		observerLib.measure("s+s-",rows,n);
	}
	if (obsOptions.find("s-s+")!=std::string::npos) {
		observerLib.measure("s-s+",rows,n);
	}
	if (obsOptions.find("ss")!=std::string::npos) {
		observerLib.measure("ss",rows,n);
	}

	//}
	return observerLib.endOfData();
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
              const std::string& targetting,
              ConcurrencyType& concurrency,
	      PsimagLite::InputValidator& io,
	      const DmrgSolverParametersType& params,
              const std::string& obsOptions)
{
	typedef Operator<RealType,MySparseMatrix> OperatorType;
	typedef Basis<RealType,MySparseMatrix> BasisType;
	typedef OperatorsBase<OperatorType,BasisType> OperatorsType;
	typedef typename OperatorType::SparseMatrixType SparseMatrixType;
	typedef BasisWithOperators<OperatorsType,ConcurrencyType> BasisWithOperatorsType; 
	typedef LeftRightSuper<BasisWithOperatorsType,BasisType> LeftRightSuperType;
	typedef ModelHelperTemplate<LeftRightSuperType,ConcurrencyType> ModelHelperType;
	typedef ModelFactory<ModelHelperType,MySparseMatrix,GeometryType,PsimagLite::PTHREADS_NAME> ModelType;
	typedef TargettingTemplate<PsimagLite::LanczosSolver,
	                           InternalProductOnTheFly,
	                           WaveFunctionTransfFactory,
	                           ModelType,
	                           ConcurrencyType,
	                           IoInputType,
	                           VectorWithOffsetTemplate> TargettingType;

	typedef DmrgSolver<InternalProductOnTheFly,TargettingType> SolverType;
	
	typedef typename TargettingType::VectorWithOffsetType VectorWithOffsetType;
	
	ModelType model(params,io,geometry,concurrency);

	 //! Read TimeEvolution if applicable:
	typedef typename TargettingType::TargettingParamsType TargettingParamsType;
	TargettingParamsType tsp(io,model);
	
	bool moreData = true;
	const std::string& datafile = params.filename;
	IoInputType dataIo(datafile);
	bool hasTimeEvolution = (targetting == "TimeStepTargetting") ? true : false;
	while (moreData) {
		try {
			moreData = !observeOneFullSweep<VectorWithOffsetType,ModelType,
			            SparseMatrixType,OperatorType,TargettingType>
			(dataIo,geometry,model,obsOptions,hasTimeEvolution,concurrency);
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
	
	ConcurrencyType concurrency(argc,argv);
	
	if (argc<2) {
		std::cerr<<"At least one argument needed\n";
		return 1;
	}
	std::string filename="";
	std::string options = "";
	if (argv[1][0]!='-') {
		std::cerr<<"WARNING: This use of the command line is deprecated.\n";
		usage(argv[0]);
		filename=argv[1];
		if (argc>2) options = argv[2];
	} else {
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
				break;
			}
		}
	}

	//Setup the Geometry
	std::vector<LabelWithKnownSizeType> labelsWithKnownSize;
	std::vector<LabelWithKnownSizeType::ValueType> vec(2,LabelWithKnownSizeType::ValueType(2,0));
	labelsWithKnownSize.push_back(LabelWithKnownSizeType("JMVALUES",vec));
	vec.resize(0);
	labelsWithKnownSize.push_back(LabelWithKnownSizeType("RAWMATRIX",vec));
	labelsWithKnownSize.push_back(LabelWithKnownSizeType("Connectors",vec));
	vec.resize(2);
	vec[0] =LabelWithKnownSizeType::ValueType(3,0);
	vec[1] =LabelWithKnownSizeType::ValueType(0,1);
	labelsWithKnownSize.push_back(LabelWithKnownSizeType( "MagneticField",vec));
	PsimagLite::InputValidator io(filename,labelsWithKnownSize);
//	IoInputType io(filename);
	GeometryType geometry(io);

	//! Read the parameters for this run
	//ParametersModelType mp(io);
	DmrgSolverParametersType dmrgSolverParams(io);
//	io.rewind();

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
			mainLoop<ModelHelperSu2,VectorWithOffsets,GroundStateTargetting,MySparseMatrixReal>
			(geometry,targetting,concurrency,io,dmrgSolverParams,options);
		} else {
			mainLoop<ModelHelperSu2,VectorWithOffset,GroundStateTargetting,MySparseMatrixReal>
			(geometry,targetting,concurrency,io,dmrgSolverParams,options);
		}
		return 0;
	}
	if (targetting=="TimeStepTargetting") { 
		mainLoop<ModelHelperLocal,VectorWithOffsets,TimeStepTargetting,MySparseMatrixComplex>
		(geometry,targetting,concurrency,io,dmrgSolverParams,options);
		return 0;
	}
	if (targetting=="DynamicTargetting") {
		mainLoop<ModelHelperLocal,VectorWithOffsets,DynamicTargetting,MySparseMatrixReal>
		(geometry,targetting,concurrency,io,dmrgSolverParams,options);
		return 0;
	}
	if (targetting=="AdaptiveDynamicTargetting") {
		mainLoop<ModelHelperLocal,VectorWithOffsets,AdaptiveDynamicTargetting,MySparseMatrixReal>
		(geometry,targetting,concurrency,io,dmrgSolverParams,options);
		return 0;
	}
	if (targetting=="CorrectionTargetting") {
		mainLoop<ModelHelperLocal,VectorWithOffsets,CorrectionTargetting,MySparseMatrixReal>
		(geometry,targetting,concurrency,io,dmrgSolverParams,options);
		return 0;
	}
	if (targetting=="MettsTargetting") { // experimental, do not use
		mainLoop<ModelHelperLocal,VectorWithOffset,MettsTargetting,MySparseMatrixReal>
		(geometry,targetting,concurrency,io,dmrgSolverParams,options);
		return 0;
	}
	mainLoop<ModelHelperLocal,VectorWithOffset,GroundStateTargetting,MySparseMatrixReal>
	(geometry,targetting,concurrency,io,dmrgSolverParams,options);
} // main

