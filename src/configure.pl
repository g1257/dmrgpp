#!/usr/bin/perl 
=pod
// BEGIN LICENSE BLOCK
Copyright (c) 2009 , UT-Battelle, LLC
All rights reserved

[DMRG++, Version 2.0.0]

*********************************************************
THE SOFTWARE IS SUPPLIED BY THE COPYRIGHT HOLDERS AND
CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED
WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A
PARTICULAR PURPOSE ARE DISCLAIMED. 

Please see full open source license included in file LICENSE.
*********************************************************


// END LICENSE BLOCK
=cut
use warnings;
use strict;

my $hasGsl = "no"; # say "no" here to remove GSL dependence

my $mpi=0;
my $platform="linux";
my $lapack="-llapack";
my $PsimagLite="../../PsimagLite/src";
my $connectors;
my $connectorValue=1.0;
my $hubbardUvalue=1.0;
my $potentialVvalue=0.0;
my ($infiniteKeptStates,$finiteLoops,$hasLoops);
my ($model,$linSize,$modelLocation,$modelFullName);
my ($geometryArgs);
my ($electrons,$momentumJ,$su2Symmetry);
my ($pthreads,$pthreadsLib)=(0,"");
my $brand= "v2.0";
my ($connectorsArgs,$connectorsArgs2,$dof,$connectors2,$connectorValue2);

my $gslLibs = " -lgsl  -lgslcblas ";
$gslLibs =" " if ($hasGsl=~/n/i);

guessPlatform();

welcome();

askQuestions();

createMakefile();

createDriver();

createObserverDriver();


sub welcome
{
	print "This is DMRG++ $brand\n";
	print "This script will help you create a Makefile, main and observer driver for your program\n";
	print "Press ENTER to continue...";
	$_=<STDIN>;
	print "\n";
	
}

sub askQuestions
{
	
	if ($platform=~/Darwin/i) {
		$lapack = " -framework Accelerate ";
	} else { # I'll assume it's linux
		$lapack = " /usr/lib64/liblapack.so.3 /usr/lib64/libblas.so.3"; 
	}

	print "Where is your PsimagLite distribution?\n";
	print "(can be obtained from https://github.com/g1257/PsimagLite)\n";
	print "Available: any\n";
	print "Default is: $PsimagLite (press ENTER): ";
	$_=<STDIN>;
	chomp;
	if ($_ eq "" or $_ eq "\n") {
		$_=$PsimagLite;
	}
	$PsimagLite=$_;

	print "Do you want to compile with MPI enabled?\n";
	print "Available: y or n\n";
	print "Default is: n (press ENTER): ";
	$_=<STDIN>;
	chomp;
	if ($_ eq "" or $_ eq "\n") {
		$_="n";
	}
	
	$mpi=1	if ($_=~/^y/i);
	
	$pthreads=0;
	$pthreadsLib="";
	if ($mpi!=1 && !($platform=~/Darwin/i)) { # cannot use both mpi and pthreads
		print "Do you want to compile with pthreads enabled?\n";
		print "Available: y or n\n";
		print "Default is: y (press ENTER): ";
		$_=<STDIN>;
		chomp;
		if ($_ eq "" or $_ eq "\n") {
			$_="y";
		}
		$pthreadsLib=" -lpthread ";
		$pthreads=1	if ($_=~/^y/i);
	}
	
	print "What model do you want to compile?\n";
	print "Available: Hubbard, Heisenberg or FeBasedSc or FeAsBasedScExtended ";
	print " or ExtendedHubbard1Orbital\n";
	print "Default is: Hubbard (press ENTER): ";
	$_=<STDIN>;
	chomp;
	if ($_ eq "" or $_ eq "\n") {
		$_="hubbard";
	}
	$model=$_;
	$modelFullName="HubbardOneBand";
	$modelFullName="HeisenbergSpinOneHalf" if ($model=~/Heisenberg/i); 
	$modelFullName="FeAsModel" if ($model=~/febasedsc/i);
	$modelFullName="TjOneOrbital" if ($model=~/tjoneorbital/i);
	$modelFullName="FeAsBasedScExtended" if ($model=~/feasbasedscextended/i);
	$modelFullName="ExtendedHubbard1Orb" if ($model=~/extendedhubbard1orb/i);
	
	$modelLocation="-IModels/$model";
	$modelLocation=" -IModels/HubbardOneBand" if ($model=~/hubbard/i);
	$modelLocation.=" -IModels/HeisenbergSpinOneHalf" if ($model=~/Heisenberg/i or $model=~/feasbasedscextended/i); 
	$modelLocation.=" -IModels/FeAsModel" if ($model=~/febasedsc/i || $model=~/feasbasedscextended/i);
	$modelLocation.=" -IModels/FeAsBasedScExtended" if ($model=~/feasbasedscextended/i);
	$modelLocation=" -IModels/TjOneOrbital" if ($model=~/tjoneorbital/i);
	$modelLocation.=" -IModels/ExtendedHubbard1Orb" if ($model=~/extendedhubbard1orb/i);
	
	if ($model=~/hubbard/i or $model=~/febasedsc/i or $model=~/tjoneorbital/i
		or $model=~/extendedhubbard1orb/i or $model=~/feasbasedscextended/i) {
		$connectors="hoppings";
	} elsif ($model=~/heisenberg/i) {
		$connectors="jvalues";
	}
	
	$connectors2="jvalues" if ($model=~/tjoneorbital/i or $model=~/feasbasedscextended/i);
	$connectors2="ninjConnectors" if ($model=~/extendedhubbard1orb/i);
	
	print "Please enter the linker flags, LDFLAGS\n";
	print "Available: Any\n";
	print "Default is: $lapack (press ENTER): ";
	$_=<STDIN>;
	chomp;
	if ($_ eq "" or $_ eq "\n") {
		$_=" $lapack ";
	}
	$lapack = $_;
}


sub createMakefile
{
	system("cp Makefile Makefile.bak") if (-r "Makefile");
	my $compiler = compilerName();
	open(FOUT,">Makefile") or die "Cannot open Makefile for writing: $!\n";
	my $usePthreadsOrNot = " ";
	$usePthreadsOrNot = " -DUSE_PTHREADS " if ($pthreads);

print FOUT<<EOF;
# DO NOT EDIT!!! Changes will be lost. Modify configure.pl instead
# This Makefile was written by configure.pl
# DMRG++ ($brand) by G.A.
# Platform: $platform
# MPI: $mpi

LDFLAGS =    $lapack  $gslLibs $pthreadsLib
CPPFLAGS = -Werror -Wall  -IEngine $modelLocation -IGeometries -I$PsimagLite $usePthreadsOrNot
EOF
if ($mpi) {
	print FOUT "CXX = mpicxx -O3 -DNDEBUG \n";
} else {
	print FOUT "CXX = $compiler  -O3 -DNDEBUG\n";
	print FOUT "#Comment out line below for debugging: \n";
	print FOUT "#CXX = $compiler -g3 \n";
}
print FOUT<<EOF;
EXENAME = dmrg
all: \$(EXENAME)

dmrg.cpp: configure.pl
	perl configure.pl

dmrg:  dmrg.o 
	\$(CXX) -o dmrg dmrg.o \$(LDFLAGS)  

correctionVectorMulti: correctionVectorMulti.o
	\$(CXX) -o correctionVectorMulti correctionVectorMulti.o \$(LDFLAGS)

# dependencies brought about by Makefile.dep
%.o: %.cpp Makefile
	\$(CXX) \$(CPPFLAGS) -c \$< 

Makefile.dep: dmrg.cpp
	\$(CXX) \$(CPPFLAGS) -MM dmrg.cpp  > Makefile.dep

observe:  observe.o
	\$(CXX) -o observe observe.o \$(LDFLAGS)

# dependencies brought about by MakefileObserver.dep
observe.o:
	\$(CXX) \$(CPPFLAGS) -c observe.cpp

MakefileObserver.dep: observe.cpp
	\$(CXX) \$(CPPFLAGS) -MM observe.cpp  > MakefileObserver.dep

clean:
	rm -f core* \$(EXENAME) *.o

include Makefile.dep
include MakefileObserver.dep
######## End of Makefile ########

EOF
	close(FOUT);
	print STDERR "File Makefile has been written\n";
}


sub computeBackwardMovements
{
	my ($directory) = @_;
	my $ret = "";
	while ($directory =~ s/\///) {
		$ret = $ret."../";
	}
	return $ret;	
}

sub createDriver
{

	system("cp dmrg.cpp dmrg.bak") if (-r "dmrg.cpp");
	open(FOUT,">dmrg.cpp") or die "Cannot open file dmrg.cpp for writing: $!\n";
	my $license=getLicense();
	my $concurrencyName = getConcurrencyName();
	my $parametersName = getParametersName();
	my $modelName = getModelName();
	my $operatorsName = getOperatorsName();
	
print FOUT<<EOF;
/* DO NOT EDIT!!! Changes will be lost. Modify configure.pl instead
 * This driver program was written by configure.pl
 * DMRG++ ($brand) by G.A.*/
#include "CrsMatrix.h"
#include "LanczosSolver.h"
#include "ChebyshevSolver.h"
#include "BlockMatrix.h"
#include "DmrgSolver.h"
#include "IoSimple.h"
#include "Operator.h"
#include "$concurrencyName.h"
#include "$modelName.h"
#include "$operatorsName.h"
#include "Geometry.h"
#ifdef USE_PTHREADS
#include "Pthreads.h"
#define PTHREADS_NAME Pthreads
#else
#include "NoPthreads.h"
#define PTHREADS_NAME NoPthreads
#endif
#include "ReflectionSymmetryEmpty.h"
#include "ModelHelperLocal.h"
#include "ModelHelperSu2.h"
#include "InternalProductOnTheFly.h"
#include "InternalProductStored.h"
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

typedef double MatrixElementType;
typedef std::complex<MatrixElementType> ComplexType;
typedef  PsimagLite::CrsMatrix<ComplexType> MySparseMatrixComplex;
typedef  PsimagLite::CrsMatrix<MatrixElementType> MySparseMatrixReal;

using namespace Dmrg;

typedef PsimagLite::$concurrencyName<MatrixElementType> ConcurrencyType;
typedef $parametersName<MatrixElementType> ParametersModelType;
typedef Geometry<MatrixElementType> GeometryType;;
typedef PsimagLite::IoSimple IoType;
typedef IoType::In IoInputType;
typedef ParametersDmrgSolver<MatrixElementType> ParametersDmrgSolverType;


template<typename ModelType,
         template<typename,typename> class InternalProductTemplate,
         typename TargettingType,
         typename MySparseMatrix>
void mainLoop3(ParametersModelType& mp,
              GeometryType& geometry,
              ParametersDmrgSolverType& dmrgSolverParams,
              ConcurrencyType& concurrency,
              IoInputType& io)
{

	//! Setup the Model
	ModelType model(mp,geometry);

	//! Read TimeEvolution if applicable:
	typedef typename TargettingType::TargettingParamsType TargettingParamsType;
	TargettingParamsType tsp(io,model);

	//! Setup the dmrg solver:
	typedef DmrgSolver<InternalProductTemplate,TargettingType> SolverType;
	SolverType dmrgSolver(dmrgSolverParams,model,concurrency,tsp);

	//! Calculate observables:
	dmrgSolver.main(geometry);
}

template<template<typename,typename,typename> class ModelHelperTemplate,
         template<typename,typename> class InternalProductTemplate,
         template<typename> class VectorWithOffsetTemplate,
         template<template<typename,typename,typename> class,
                  template<typename,typename> class,
                  template<typename,typename> class,
                  typename,typename,typename,
                  template<typename> class> class TargettingTemplate,
         typename MySparseMatrix>
void mainLoop2(ParametersModelType& mp,
              GeometryType& geometry,
              ParametersDmrgSolverType& dmrgSolverParams,
              ConcurrencyType& concurrency,
              IoInputType& io)
{
	typedef ReflectionSymmetryEmpty<MatrixElementType,MySparseMatrix> ReflectionSymmetryType;
	typedef Operator<MatrixElementType,MySparseMatrix> OperatorType;
	typedef Basis<MatrixElementType,MySparseMatrix> BasisType;
	typedef $operatorsName<OperatorType,BasisType> OperatorsType;
	typedef BasisWithOperators<OperatorsType,ConcurrencyType> BasisWithOperatorsType;
	typedef LeftRightSuper<BasisWithOperatorsType,BasisType> LeftRightSuperType;
	typedef ModelHelperTemplate<LeftRightSuperType,ReflectionSymmetryType,ConcurrencyType> ModelHelperType;
	typedef $modelName<ModelHelperType,MySparseMatrix,GeometryType,PsimagLite::PTHREADS_NAME> ModelType;
	
	if (dmrgSolverParams.options.find("ChebyshevSolver")!=std::string::npos) {
			typedef TargettingTemplate<PsimagLite::ChebyshevSolver,
			                           InternalProductTemplate,
			                           WaveFunctionTransfFactory,
			                           ModelType,
			                           ConcurrencyType,
			                           IoType,
			                           VectorWithOffsetTemplate
			                           > TargettingType;
			mainLoop3<ModelType,InternalProductTemplate,TargettingType,MySparseMatrix>
			(mp,geometry,dmrgSolverParams,concurrency,io);
	} else {
		typedef TargettingTemplate<PsimagLite::LanczosSolver,
			                           InternalProductTemplate,
			                           WaveFunctionTransfFactory,
			                           ModelType,
			                           ConcurrencyType,
			                           IoType,
			                           VectorWithOffsetTemplate
			                           > TargettingType;
			mainLoop3<ModelType,InternalProductTemplate,TargettingType,MySparseMatrix>
			(mp,geometry,dmrgSolverParams,concurrency,io);
	}
}

template<template<typename,typename,typename> class ModelHelperTemplate,
         template<typename> class VectorWithOffsetTemplate,
         template<template<typename,typename,typename> class,
                  template<typename,typename> class,
                  template<typename,typename> class,
                  typename,typename,typename,
                  template<typename> class> class TargettingTemplate,
         typename MySparseMatrix>
void mainLoop(ParametersModelType& mp,
              GeometryType& geometry,
              ParametersDmrgSolverType& dmrgSolverParams,
              ConcurrencyType& concurrency,
              IoInputType& io)
{
	if (dmrgSolverParams.options.find("InternalProductStored")!=std::string::npos) {
		mainLoop2<ModelHelperTemplate,
		         InternalProductStored,
		         VectorWithOffsetTemplate,
		         TargettingTemplate,
		         MySparseMatrix>(mp,geometry,dmrgSolverParams,concurrency,io);
	} else {
 		mainLoop2<ModelHelperTemplate,
		         InternalProductOnTheFly,
		         VectorWithOffsetTemplate,
		         TargettingTemplate,
		         MySparseMatrix>(mp,geometry,dmrgSolverParams,concurrency,io);
	}
}

int main(int argc,char *argv[])
{
	//! setup distributed parallelization
	ConcurrencyType concurrency(argc,argv);

	// print license
	std::string license = $license;
	if (concurrency.root()) std::cerr<<license;

	//Setup the Geometry
	IoInputType io(argv[1]);
	GeometryType geometry(io);

	//! Read the parameters for this run
	ParametersModelType mp(io);
	ParametersDmrgSolver<MatrixElementType> dmrgSolverParams(io);

	bool su2=false;
	if (dmrgSolverParams.options.find("useSu2Symmetry")!=std::string::npos) su2=true;
	std::string targetting="GroundStateTargetting";
	const char *targets[]={"TimeStepTargetting","DynamicTargetting","AdaptiveDynamicTargetting",
                     "CorrectionVectorTargetting","CorrectionTargetting","MettsTargetting"};
	size_t totalTargets = 6;
	for (size_t i = 0;i<totalTargets;++i)
		if (dmrgSolverParams.options.find(targets[i])!=std::string::npos) targetting=targets[i];

	if (targetting!="GroundStateTargetting" && su2) throw std::runtime_error("SU(2)"
 		" supports only GroundStateTargetting for now (sorry!)\\n");
	if (su2) {
		if (dmrgSolverParams.targetQuantumNumbers[2]>0) { 
			mainLoop<ModelHelperSu2,VectorWithOffsets,GroundStateTargetting,
				MySparseMatrixReal>(mp,geometry,dmrgSolverParams,concurrency,io);
		} else {
			mainLoop<ModelHelperSu2,VectorWithOffset,GroundStateTargetting,
				MySparseMatrixReal>(mp,geometry,dmrgSolverParams,concurrency,io);
		}
		return 0;
	}
	if (targetting=="TimeStepTargetting") { 
		mainLoop<ModelHelperLocal,VectorWithOffsets,TimeStepTargetting,
			MySparseMatrixComplex>(mp,geometry,dmrgSolverParams,concurrency,io);
			return 0;
	}
	if (targetting=="DynamicTargetting") {
		mainLoop<ModelHelperLocal,VectorWithOffsets,DynamicTargetting,
			MySparseMatrixReal>(mp,geometry,dmrgSolverParams,concurrency,io);
			return 0;
	}
	if (targetting=="AdaptiveDynamicTargetting") {
		mainLoop<ModelHelperLocal,VectorWithOffsets,AdaptiveDynamicTargetting,
			MySparseMatrixReal>(mp,geometry,dmrgSolverParams,concurrency,io);
			return 0;
	}
	if (targetting=="CorrectionVectorTargetting") {
		mainLoop<ModelHelperLocal,VectorWithOffsets,CorrectionVectorTargetting,
			MySparseMatrixReal>(mp,geometry,dmrgSolverParams,concurrency,io);
			return 0;
	}
	if (targetting=="CorrectionTargetting") {
		mainLoop<ModelHelperLocal,VectorWithOffsets,CorrectionTargetting,
			MySparseMatrixReal>(mp,geometry,dmrgSolverParams,concurrency,io);
			return 0;
	}
	if (targetting=="MettsTargetting") {
		mainLoop<ModelHelperLocal,VectorWithOffset,MettsTargetting,
			MySparseMatrixReal>(mp,geometry,dmrgSolverParams,concurrency,io);
			return 0;
	}
	mainLoop<ModelHelperLocal,VectorWithOffset,GroundStateTargetting,
		MySparseMatrixReal>(mp,geometry,dmrgSolverParams,concurrency,io);
}

EOF
	close(FOUT);
	print STDERR "File dmrg.cpp has been written\n";

}

sub getLicense
{
	open(THISFILE,"$0") or return " ";
	while(<THISFILE>) {
		last if (/BEGIN LICENSE BLOCK/);
	}
	my $l="";
	while(<THISFILE>) {
		last if (/END LICENSE BLOCK/);
		chomp;
		s/\"/\\\"/g;
		$l = "$l\"$_\\n\"\n";
	}
	close(THISFILE);
	return $l;
}

sub createObserverDriver
{
	my $observerDriver="observe.cpp";
	my $license=getLicense();
	my $concurrencyName = getConcurrencyName();
	my $parametersName = getParametersName();
	my $modelName = getModelName();
	my $operatorsName = getOperatorsName();

	system("cp observe.cpp observe.bak") if (-e "observe.cpp");	
	open(OBSOUT,">$observerDriver") or die "Cannot open file $observerDriver for writing: $!\n";
print OBSOUT<<EOF;
/* DO NOT EDIT!!! Changes will be lost. Modify configure.pl instead
 * This driver program was written by configure.pl
 * DMRG++ ($brand) by G.A.*/
	
#include "Observer.h"
#include "ObservableLibrary.h"
#include "IoSimple.h"
#include "$modelName.h" 
#include "$operatorsName.h" 
#include "$concurrencyName.h" 
#include "Geometry.h" 
#include "CrsMatrix.h"
#include "ReflectionSymmetryEmpty.h"
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

using namespace Dmrg;

typedef double RealType;
typedef std::complex<RealType> ComplexType;

typedef  PsimagLite::CrsMatrix<ComplexType> MySparseMatrixComplex;
typedef  PsimagLite::CrsMatrix<RealType> MySparseMatrixReal;
typedef Geometry<RealType> GeometryType;
typedef  $parametersName<RealType> ParametersModelType; 
typedef PsimagLite::$concurrencyName<RealType> ConcurrencyType;
typedef PsimagLite::IoSimple::In IoInputType;

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
	std::string sSweeps = "sweeps=";
	std::string::size_type begin = obsOptions.find(sSweeps);
	if (begin != std::string::npos) {
		std::string sTmp = obsOptions.substr(begin+sSweeps.length(),std::string::npos);
		//std::cout<<"sTmp="<<sTmp<<"\\n";
		n = atoi(sTmp.c_str());
	}
	ObservableLibraryType observerLib(io,n,hasTimeEvolution,model,concurrency,verbose);
	
	bool ot = false;
	if (obsOptions.find("ot") || obsOptions.find("time")) ot = true;
	if (hasTimeEvolution && ot) {
		observerLib.measureTime("superDensity");
		observerLib.measureTime("nupNdown");
		observerLib.measureTime("nup+ndown");
	}

	if (hasTimeEvolution) observerLib.setBrackets("time","time");
EOF
	if  ($modelName=~/heisenberg/i) {
	} else {
print OBSOUT<<EOF; 

	if (obsOptions.find("cc")!=std::string::npos) {
		observerLib.measure("cc",n/2,n);
	}
	
	if (obsOptions.find("nn")!=std::string::npos) {
		observerLib.measure("nn",n/2,n);
	}
EOF
	}
print OBSOUT<<EOF;
	if (obsOptions.find("szsz")!=std::string::npos) {
		observerLib.measure("szsz",n/2,n);
	}
EOF
	if  ($modelName=~/febasedsc/i or $modelName=~/feasbasedscextended/i) {
print print OBSOUT<<EOF;
	if (obsOptions.find("dd")!=std::string::npos && 
			geometry.label(0).find("ladder")==std::string::npos) {
		observerLib.measure("dd",n/2,n);
	}
	

	// FOUR-POINT DELTA-DELTA^DAGGER:
	if (obsOptions.find("dd4")!=std::string::npos &&
		geometry.label(0).find("ladder")!=std::string::npos) {
		observerLib.measure("dd4",n/2,n);
	} // if dd4
EOF
	}
	if  ($modelName=~/heisenberg/i) {
	print print OBSOUT<<EOF;
		if (obsOptions.find("s+s-")!=std::string::npos) {
			observerLib.measure("s+s-",n/2,n);
		}
		if (obsOptions.find("s-s+")!=std::string::npos) {
			observerLib.measure("s-s+",n/2,n);
		}
		if (obsOptions.find("ss")!=std::string::npos) {
			observerLib.measure("ss",n/2,n);
		}
EOF
	}
	print OBSOUT<<EOF;
	return observerLib.endOfData();
}

template<template<typename,typename,typename> class ModelHelperTemplate,
         template<typename> class VectorWithOffsetTemplate,
         template<template<typename,typename,typename> class,
                  template<typename,typename> class,
                  template<typename,typename> class,
                  typename,typename,typename,
         template<typename> class> class TargettingTemplate,
         typename MySparseMatrix>
void mainLoop(ParametersModelType& mp,
              GeometryType& geometry,
              const std::string& targetting,
              ConcurrencyType& concurrency,
              IoInputType& io,
              const std::string& datafile,
              const std::string& obsOptions)
{
	typedef ReflectionSymmetryEmpty<RealType,MySparseMatrix> ReflectionSymmetryType;
	typedef Operator<RealType,MySparseMatrix> OperatorType;
	typedef Basis<RealType,MySparseMatrix> BasisType;
	typedef $operatorsName<OperatorType,BasisType> OperatorsType;
	typedef typename OperatorType::SparseMatrixType SparseMatrixType;
	typedef BasisWithOperators<OperatorsType,ConcurrencyType> BasisWithOperatorsType; 
	typedef LeftRightSuper<BasisWithOperatorsType,BasisType> LeftRightSuperType;
	typedef ModelHelperTemplate<LeftRightSuperType,ReflectionSymmetryType,ConcurrencyType> ModelHelperType;
	typedef $modelName<ModelHelperType,MySparseMatrix,GeometryType,PsimagLite::PTHREADS_NAME> ModelType;
	typedef TargettingTemplate<PsimagLite::LanczosSolver,
	                           InternalProductOnTheFly,
	                           WaveFunctionTransfFactory,
	                           ModelType,
	                           ConcurrencyType,
	                           IoInputType,
	                           VectorWithOffsetTemplate> TargettingType;

	typedef DmrgSolver<InternalProductOnTheFly,TargettingType> SolverType;
	
	typedef typename TargettingType::VectorWithOffsetType VectorWithOffsetType;
	
	ModelType model(mp,geometry);

	 //! Read TimeEvolution if applicable:
	typedef typename TargettingType::TargettingParamsType TargettingParamsType;
	TargettingParamsType tsp(io,model);
	
	bool moreData = true;
	IoInputType dataIo(datafile);
	bool hasTimeEvolution = (targetting == "TimeStepTargetting") ? true : false;
	while (moreData) {
		try {
			moreData = !observeOneFullSweep<VectorWithOffsetType,ModelType,
			            SparseMatrixType,OperatorType,TargettingType>
			(dataIo,geometry,model,obsOptions,hasTimeEvolution,concurrency);
		} catch (std::exception& e) {
			std::cerr<<"CAUGHT: "<<e.what();
			std::cerr<<"There's no more data\\n";
			break;
		}

		if (!hasTimeEvolution) break;
	}
}

int main(int argc,char *argv[])
{
	
	using namespace Dmrg;
	
	typedef   PsimagLite::$concurrencyName<RealType> ConcurrencyType; 
	
	ConcurrencyType concurrency(argc,argv);
	
	if (argc<2) {
		std::string s = "The observer driver takes at least one  argumentw: \\n";
		s = s + "(i) the name of the input file";
		s = s+ " and, optionally, (ii) a comma-separated list of options of what to compute\\n";
		throw std::runtime_error(s);
	}
	std::string options = "";
	if (argc>2) options = argv[2];
	
	//Setup the Geometry
	IoInputType io(argv[1]);
	GeometryType geometry(io);

	//! Read the parameters for this run
	ParametersModelType mp(io);
	ParametersDmrgSolver<RealType> dmrgSolverParams(io);

	bool su2=false;
	if (dmrgSolverParams.options.find("useSu2Symmetry")!=std::string::npos) su2=true;
	std::string targetting="GroundStateTargetting";
	const char *targets[]={"TimeStepTargetting","DynamicTargetting","AdaptiveDynamicTargetting",
                     "CorrectionVectorTargetting","CorrectionTargetting","MettsTargetting"};
	size_t totalTargets = 6;
	for (size_t i = 0;i<totalTargets;++i)
		if (dmrgSolverParams.options.find(targets[i])!=std::string::npos) targetting=targets[i];
	if (targetting!="GroundStateTargetting" && su2) throw std::runtime_error("SU(2)"
 		" supports only GroundStateTargetting for now (sorry!)\\n");
	
	if (su2) {
		if (dmrgSolverParams.targetQuantumNumbers[2]>0) { 
			mainLoop<ModelHelperSu2,VectorWithOffsets,GroundStateTargetting,MySparseMatrixReal>
			(mp,geometry,targetting,concurrency,io,dmrgSolverParams.filename,options);
		} else {
			mainLoop<ModelHelperSu2,VectorWithOffset,GroundStateTargetting,MySparseMatrixReal>
			(mp,geometry,targetting,concurrency,io,dmrgSolverParams.filename,options);
		}
		return 0;
	}
	if (targetting=="TimeStepTargetting") { 
		mainLoop<ModelHelperLocal,VectorWithOffsets,TimeStepTargetting,MySparseMatrixComplex>
		(mp,geometry,targetting,concurrency,io,dmrgSolverParams.filename,options);
		return 0;
	}
	if (targetting=="DynamicTargetting") {
		mainLoop<ModelHelperLocal,VectorWithOffsets,DynamicTargetting,MySparseMatrixReal>
		(mp,geometry,targetting,concurrency,io,dmrgSolverParams.filename,options);
		return 0;
	}
	if (targetting=="AdaptiveDynamicTargetting") {
		mainLoop<ModelHelperLocal,VectorWithOffsets,AdaptiveDynamicTargetting,MySparseMatrixReal>
		(mp,geometry,targetting,concurrency,io,dmrgSolverParams.filename,options);
		return 0;
	}
	if (targetting=="CorrectionTargetting") {
		mainLoop<ModelHelperLocal,VectorWithOffsets,CorrectionTargetting,MySparseMatrixReal>
		(mp,geometry,targetting,concurrency,io,dmrgSolverParams.filename,options);
		return 0;
	}
	if (targetting=="MettsTargetting") { // experimental, do not use
		mainLoop<ModelHelperLocal,VectorWithOffset,MettsTargetting,MySparseMatrixReal>
		(mp,geometry,targetting,concurrency,io,dmrgSolverParams.filename,options);
		return 0;
	}
	mainLoop<ModelHelperLocal,VectorWithOffset,GroundStateTargetting,MySparseMatrixReal>
	(mp,geometry,targetting,concurrency,io,dmrgSolverParams.filename,options);
} // main

EOF
	close(OBSOUT);
	print STDERR "File observe.cpp has been written\n";
}

sub guessPlatform
{
	$platform="Linux";
	$platform="Darwin" if (isAMac());
}

sub isAMac
{
	open(PIPE,"uname -a |grep -i Darwin") or return 0;
	$_=<PIPE>;
	close(PIPE);
	
	return 1 unless ($_ eq "" or $_ eq "\n");
	return 0;
}


sub getParametersName
{
	my $parametersName="Parameters$model";
	if ($model=~/hubbard/i) {
		$parametersName = "ParametersModelHubbard";
	} elsif ($model=~/heisenberg/i) {
		$parametersName = "ParametersModelHeisenberg";
	} elsif ($model=~/febasedsc/i) {
		$parametersName = "ParametersModelFeAs";
	} elsif ($model=~/feasbasedscextended/i) {
		$parametersName = "ParametersModelFeAs";
	} elsif ($model=~/tjoneorbital/i) {
		$parametersName = "ParametersTjOneOrbital";
	}
	return $parametersName;
}

sub getConcurrencyName()
{
	my $concurrencyName = "UNKNOWN";

	if ($mpi) {
		$concurrencyName="ConcurrencyMpi";
	} else {
		$concurrencyName="ConcurrencySerial";
	}
	return $concurrencyName;
}

sub getModelName()
{
	my $modelName = $model;

	if ($model=~/heisenberg/i) {
		$modelName="ModelHeisenberg";
	} elsif ($model=~/febasedsc/i) {
		$modelName = "ModelFeBasedSc";
	} elsif ($model=~/feasbasedscextended/i) {
		$modelName = "FeAsBasedScExtended";
	} elsif ($model=~/tjoneorbital/i) {
		$modelName = "TjOneOrbital";
	} elsif ($model=~/extendedhubbard1orb/i) {
		$modelName = "ExtendedHubbard1Orb";
	} elsif ($model=~/hubbard/i) {
		$modelName = "ModelHubbard"; # after extended
	}
	return $modelName;
}

sub getOperatorsName()
{
	my $operatorsName = "Operators$model";

	if ($model=~/extendedhubbard1orb/i) {
		$operatorsName = "OpsExtendedHubbard1Orb";
	} elsif ($model=~/hubbard/i) {
		$operatorsName = "OperatorsHubbard";
	} elsif ($model=~/heisenberg/i) {
		$operatorsName = "OperatorsHeisenberg";
	} elsif ($model=~/febasedsc/i) {
		$operatorsName="OperatorsFeAs";
	} elsif ($model=~/feasbasedscextended/i) {
		$operatorsName="OperatorsFeAsExtended";
	} elsif ($model=~/tjoneorbital/i) {
		$operatorsName="OperatorsTjOneOrbital";
	}
	return $operatorsName;
}

sub compilerName
{
	return "g++";
	my @tryThis = ("g++","g++4");
	my $ret;
	my $compiler;
	foreach my $comp (@tryThis) {
		my $ret = system("$comp &>2 /dev/null");
		if ($ret==0) {
			$compiler = $comp;
			last;
		} 
			
	}
	return $compiler if defined $compiler; 
	die "$0: No suitable compiler found\n";
}

