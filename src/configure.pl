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
my $targetting;
my $DynamicTargetting = "DynamicTargetting";
#$DynamicTargetting = "DynamicTargettingEmpty" if ($hasGsl=~/n/i);

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
		$lapack = " -llapack -lblas"; 
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
	print " or TjOneOrbital or ExtendedHubbard1Orbital\n";
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
	
	$modelLocation="";
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
	
	print "What should the observer driver target?\n";
	print "Available: GroundStateTargetting TimeStepTargetting $DynamicTargetting\n";
	print "Default is: GroundStateTargetting (press ENTER): ";
	$_=<STDIN>;
	chomp;
	if ($_ eq "" or $_ eq "\n") {
		$_="GroundStateTargetting";
	}
	$targetting = $_;
	
	
}


sub createMakefile
{
	system("cp Makefile Makefile.bak") if (-r "Makefile");
	my $compiler = compilerName();
	my $headerFiles = join(' ', glob("Engine/*.h Models/*/*.h Geometries/*.h"));
	my @litProgFiles = glob("Engine/*.w Models/*/*.w Geometries/*.w");
	my $litProgTargets = getLitProgTargets(\@litProgFiles);
	open(FOUT,">Makefile") or die "Cannot open Makefile for writing: $!\n";
print FOUT<<EOF;
# DO NOT EDIT!!! Changes will be lost. Modify configure.pl instead
# This Makefile was written by configure.pl
# DMRG++ ($brand) by G.A.
# Platform: $platform
# MPI: $mpi

LDFLAGS =    $lapack  $gslLibs $pthreadsLib
CPPFLAGS = -Werror -Wall  -IEngine $modelLocation -IGeometries -I$PsimagLite
EOF
if ($mpi) {
	print FOUT "CXX = mpicxx -O2 -DNDEBUG \n";
} else {
	print FOUT "CXX = $compiler -pg -O2 -DNDEBUG\n";
}
print FOUT<<EOF;
all: \$(EXENAME)
HEADERSH = $headerFiles

all: dmrg

dmrg:  \$(HEADERSH) 
	\$(CXX) -o dmrg  \$(CPPFLAGS)  dmrg.cpp \$(LDFLAGS)  

observe:  \$(HEADERSH)
	\$(CXX) -o observe \$(CPPFLAGS) observe.cpp \$(LDFLAGS)

lanczos: \$(HEADERSH)
	\$(CXX) -o lanczos \$(CPPFLAGS) lanczos.cpp \$(LDFLAGS)

clean:
	rm -f core* \$(EXENAME) *.o *.ii *.tt

$litProgTargets
######## End of Makefile ########

EOF
	close(FOUT);
	print STDERR "File Makefile has been written\n";
}

sub getLitProgTargets
{
	my ($array)=@_;
	my $x = "";
	my $litProgTool = "nuweb.pl -v -l  -s  -d ";
	foreach my $f (@$array) {
		my $fh = $f;
		$fh =~ s/\.w$/\.h/;
		$x = $x."$fh: $f\n";
		my $dir = $f;
		$dir =~ s/\/[^\/]+$/\//;
		my $fnd = $f;
		$fnd =~ s/$dir//;
		my $dirChange = computeBackwardMovements($dir);
		$x = $x."\t cd $dir; $dirChange$litProgTool $fnd\n";
		$x = $x."\n";
	}
	return $x;
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
	my $pthreadsName = getPthreadsName();
	my $modelName = getModelName();
	my $operatorsName = getOperatorsName();
	
print FOUT<<EOF;
/* DO NOT EDIT!!! Changes will be lost. Modify configure.pl instead
 * This driver program was written by configure.pl
 * DMRG++ ($brand) by G.A.*/
#include "CrsMatrix.h"
#include "LanczosSolver.h"
#include "BlockMatrix.h"
#include "DmrgSolver.h"
#include "IoSimple.h"
#include "Operator.h"
#include "$concurrencyName.h"
#include "$modelName.h"
#include "$operatorsName.h"
#include "Geometry.h"
#include "$pthreadsName.h"
#include "ReflectionSymmetryEmpty.h"
#include "ModelHelperLocal.h"
#include "ModelHelperSu2.h"
#include "InternalProductOnTheFly.h"
#include "InternalProductStored.h"
#include "GroundStateTargetting.h"
#include "TimeStepTargetting.h"
#include "$DynamicTargetting.h"
#include "VectorWithOffset.h"
#include "VectorWithOffsets.h"
#include "BasisWithOperators.h"
#include "LeftRightSuper.h"

typedef double MatrixElementType;
typedef std::complex<MatrixElementType> ComplexType;
typedef  PsimagLite::CrsMatrix<ComplexType> MySparseMatrixComplex;
typedef  PsimagLite::CrsMatrix<MatrixElementType> MySparseMatrixReal;

using namespace Dmrg;

typedef double Field;
//typedef std::vector<int> MyBlock;
//typedef BlockMatrix<MatrixElementType,MySparseMatrix> MySparseBlockMatrix;
//typedef BlockMatrix<MatrixElementType,PsimagLite::Matrix<MatrixElementType> > MyDenseBlockMatrix;
typedef PsimagLite::$concurrencyName<MatrixElementType> MyConcurrency;
typedef $parametersName<MatrixElementType> ParametersModelType;
typedef Geometry<MatrixElementType> GeometryType;;
typedef  PsimagLite::IoSimple MyIo;

template<
	typename ParametersModelType,
	typename GeometryType,
	typename ParametersSolverType,
	typename ConcurrencyType,
	typename IoInputType,
	template<typename,typename,typename,template<typename> class> class ModelTemplate,
	template<typename,typename,typename> class ModelHelperTemplate,
	template<typename,typename> class InternalProductTemplate,
	template<typename> class VectorWithOffsetTemplate,
	template<template<typename,typename,typename> class,
		template<typename,typename> class,
		template<typename,typename> class,
		typename,typename,typename,
		template<typename> class> class TargettingTemplate,
	typename MySparseMatrix
>
void mainLoop(ParametersModelType& mp,GeometryType& geometry,ParametersSolverType& dmrgSolverParams,
		ConcurrencyType& concurrency, IoInputType& io,const std::string& targettingString)
{
	typedef ReflectionSymmetryEmpty<MatrixElementType,MySparseMatrix> ReflectionSymmetryType;
	typedef Operator<MatrixElementType,MySparseMatrix> OperatorType;
	typedef Basis<MatrixElementType,MySparseMatrix> BasisType;
	typedef $operatorsName<OperatorType,BasisType> OperatorsType;
	typedef BasisWithOperators<OperatorsType,ConcurrencyType> BasisWithOperatorsType;
	typedef LeftRightSuper<BasisWithOperatorsType,BasisType> LeftRightSuperType;
	typedef ModelHelperTemplate<LeftRightSuperType,ReflectionSymmetryType,MyConcurrency> ModelHelperType;
	typedef ModelTemplate<ModelHelperType,MySparseMatrix,GeometryType,PsimagLite::$pthreadsName> ModelType;
	
	typedef DmrgSolver<
			InternalProductTemplate,
			ModelHelperTemplate,
			ModelType,
			MyIo,
			TargettingTemplate,
			VectorWithOffsetTemplate
		> SolverType;
	
	//! Setup the Model
	ModelType model(mp,geometry);

	//! Read TimeEvolution if applicable:
	typedef typename SolverType::TargettingType TargettingType;
	typedef typename TargettingType::TargettingParamsType TargettingParamsType;
	TargettingParamsType tsp(io,model);

	//! Setup the dmrg solver:
	SolverType dmrgSolver(dmrgSolverParams,model,concurrency,tsp);

	//! Calculate observables:
	dmrgSolver.main(geometry);
}

int main(int argc,char *argv[])
{
	//! setup distributed parallelization
	MyConcurrency concurrency(argc,argv);
	
	//Setup the Geometry
	typedef PsimagLite::IoSimple::In IoInputType;
	IoInputType io(argv[1]);
	GeometryType geometry(io);

	//! Read the parameters for this run
	ParametersModelType mp(io);
	ParametersDmrgSolver<MatrixElementType> dmrgSolverParams(io);

	// print license
	std::string license = $license;
	if (concurrency.root()) std::cerr<<license;

	bool su2=false;
	if (dmrgSolverParams.options.find("useSu2Symmetry")!=std::string::npos) su2=true;
	std::string targetting="GroundStateTargetting";
	if (dmrgSolverParams.options.find("TimeStepTargetting")!=std::string::npos) targetting="TimeStepTargetting";
	if (dmrgSolverParams.options.find("DynamicTargetting")!=std::string::npos) targetting="DynamicTargetting";
	if (targetting!="GroundStateTargetting" && su2) throw std::runtime_error("SU(2)"
 		" supports only GroundStateTargetting for now (sorry!)\\n");
	if (su2) {
		if (dmrgSolverParams.targetQuantumNumbers[2]>0) { 
			mainLoop<ParametersModelType,GeometryType,ParametersDmrgSolver<MatrixElementType>,MyConcurrency,
				IoInputType,
				$modelName,ModelHelperSu2,InternalProductOnTheFly,VectorWithOffsets,GroundStateTargetting,
				MySparseMatrixReal>
			(mp,geometry,dmrgSolverParams,concurrency,io,targetting);
		} else {
			mainLoop<ParametersModelType,GeometryType,ParametersDmrgSolver<MatrixElementType>,MyConcurrency,
				IoInputType,
				$modelName,ModelHelperSu2,InternalProductOnTheFly,VectorWithOffset,GroundStateTargetting,
				MySparseMatrixReal>
			(mp,geometry,dmrgSolverParams,concurrency,io,targetting);
		}
		return 0;
	}
	if (targetting=="TimeStepTargetting") { 
		mainLoop<ParametersModelType,GeometryType,ParametersDmrgSolver<MatrixElementType>,MyConcurrency,
			IoInputType,
			$modelName,ModelHelperLocal,InternalProductOnTheFly,VectorWithOffsets,TimeStepTargetting,
			MySparseMatrixComplex>
			(mp,geometry,dmrgSolverParams,concurrency,io,targetting);
			return 0;
	}
	if (targetting=="DynamicTargetting") {
		mainLoop<ParametersModelType,GeometryType,ParametersDmrgSolver<MatrixElementType>,MyConcurrency,
			IoInputType,
			$modelName,ModelHelperLocal,InternalProductOnTheFly,VectorWithOffsets,DynamicTargetting,
			MySparseMatrixReal>
			(mp,geometry,dmrgSolverParams,concurrency,io,targetting);
			return 0;
	}
	
	mainLoop<ParametersModelType,GeometryType,ParametersDmrgSolver<MatrixElementType>,MyConcurrency,
		IoInputType,
		$modelName,ModelHelperLocal,InternalProductOnTheFly,VectorWithOffset,GroundStateTargetting,
		MySparseMatrixReal>
		(mp,geometry,dmrgSolverParams,concurrency,io,targetting);
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
		chomp;
		s/\"/\\\"/g;
		$l = "$l\"$_\\n\"\n";
		last if (/END LICENSE BLOCK/);
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
	my $pthreadsName = getPthreadsName();
	my $modelName = getModelName();
	my $operatorsName = getOperatorsName();
	my $chooseRealOrComplexForObservables = "typedef RealType FieldType;\n";
	if ($targetting=~/timestep/i) {
		$chooseRealOrComplexForObservables = "typedef ComplexType FieldType;\n";
	}

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
#include "$pthreadsName.h" 
#include "ModelHelperLocal.h"
#include "ModelHelperSu2.h"
#include "InternalProductOnTheFly.h"
#include "VectorWithOffset.h"
#include "VectorWithOffsets.h"
#include "GroundStateTargetting.h"
#include "DmrgSolver.h" // only used for types
#include "TimeStepTargetting.h" // only used for types
#include "GroundStateTargetting.h" // only used for types
#include "BasisWithOperators.h"
#include "LeftRightSuper.h"

using namespace Dmrg;

typedef double RealType;
typedef std::complex<RealType> ComplexType;
$chooseRealOrComplexForObservables

typedef PsimagLite::$concurrencyName<RealType> MyConcurrency;
typedef PsimagLite::IoSimple::In IoInputType;

template<typename ConcurrencyType,typename VectorWithOffsetType,typename ModelType,typename SparseMatrixType,
typename OperatorType,typename TargettingType,typename GeometryType>
bool observeOneFullSweep(IoInputType& io,
	const GeometryType& geometry,const ModelType& model,const std::string& obsOptions,
		bool hasTimeEvolution,ConcurrencyType& concurrency)
{
	bool verbose = false;
	size_t n=geometry.numberOfSites();
	typedef Observer<FieldType,VectorWithOffsetType,ModelType,IoInputType> 
		ObserverType;
	typedef ObservableLibrary<ObserverType,TargettingType> ObservableLibraryType;
	ObservableLibraryType observerLib(io,n,hasTimeEvolution,model,concurrency,verbose);
	
	my $ot = 0;
	if (obsOptions.find("ot") || obsOptions.find("time")) $ot = 1;
	if (hasTimeEvolution && $ot) {
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

	template<
	typename ParametersModelType,
	typename GeometryType,
	typename ConcurrencyType,
	typename IoInputType,
	template<typename,typename,typename,template<typename> class> class ModelTemplate,
	template<typename,typename,typename> class ModelHelperTemplate,
	template<typename,typename> class InternalProductTemplate,
	template<typename> class VectorWithOffsetTemplate,
	template<template<typename,typename,typename> class,
		template<typename,typename> class,
		template<typename,typename> class,
		typename,typename,typename,
		template<typename> class> class TargettingTemplate,
	typename MySparseMatrix
>
void mainLoop(ParametersModelType& mp,GeometryType& geometry,bool hasTimeEvolution,
		ConcurrencyType& concurrency, IoInputType& io,const std::string& datafile,
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
	typedef ModelTemplate<ModelHelperType,MySparseMatrix,GeometryType,PsimagLite::$pthreadsName> ModelType;
	
	typedef DmrgSolver<
                        InternalProductTemplate,
			ModelHelperTemplate,
                        ModelType,
                        PsimagLite::IoSimple,
                        TargettingTemplate,
                        VectorWithOffsetTemplate
                > SolverType; // only used for types
	
	typedef typename SolverType::TargettingType TargettingType;
	typedef typename TargettingType::VectorWithOffsetType VectorWithOffsetType;
	
	ModelType model(mp,geometry);

	 //! Read TimeEvolution if applicable:
	typedef typename SolverType::TargettingType TargettingType;
	typedef typename TargettingType::TargettingParamsType TargettingParamsType;
	TargettingParamsType tsp(io,model);
	
	//size_t n=geometry.numberOfSites();
	bool moreData = true;
	IoInputType dataIo(datafile);
	while (moreData) {
		try {
			moreData = !observeOneFullSweep<ConcurrencyType,VectorWithOffsetType,ModelType,
				SparseMatrixType,OperatorType,TargettingType,GeometryType>
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
	
	typedef   PsimagLite::$concurrencyName<RealType> MyConcurrency; 
	
	MyConcurrency concurrency(argc,argv);
	
	if (argc<2) {
		std::string s = "The observer driver takes at least one  argumentw: \\n";
		s = s + "(i) the name of the input file";
		s = s+ " and, optionally, (ii) a comma-separated list of options of what to compute\\n";
		throw std::runtime_error(s);
	}
	std::string options = "";
	if (argc>2) options = argv[2];
	
	//Setup the Geometry
	typedef Geometry<RealType> GeometryType;
	IoInputType io(argv[1]);
	GeometryType geometry(io);

	//! Read the parameters for this run
	typedef  $parametersName<RealType> ParametersModelType; 
	ParametersModelType mp(io);
	ParametersDmrgSolver<FieldType> dmrgSolverParams(io);

	bool hasTimeEvolution=false;
	if (dmrgSolverParams.options.find("TimeStepTargetting")!=std::string::npos) hasTimeEvolution=true;
	
	// FIXME: See if we need VectorWithOffsets sometimes
	// FIXME: Does it make sense to have ModelHelperSu2 here sometimes?
	typedef PsimagLite::CrsMatrix<RealType> MySparseMatrixReal;
	typedef PsimagLite::CrsMatrix<FieldType> MySparseMatrixComplex;
	if (hasTimeEvolution)
		mainLoop<ParametersModelType,GeometryType,MyConcurrency,IoInputType,$modelName,
			ModelHelperLocal,InternalProductOnTheFly,VectorWithOffsets,
			$targetting,MySparseMatrixComplex>(mp,geometry,hasTimeEvolution,concurrency,io,dmrgSolverParams.filename,options);
	else
		mainLoop<ParametersModelType,GeometryType,MyConcurrency,IoInputType,$modelName,
			ModelHelperLocal,InternalProductOnTheFly,VectorWithOffset,
			$targetting,MySparseMatrixReal>(mp,geometry,hasTimeEvolution,concurrency,io,dmrgSolverParams.filename,options);
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
	my $parametersName="UNKNOWN";
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

sub getPthreadsName
{
	my $pthreadsName = "UNKNOWN";
	if ($pthreads) {
		$pthreadsName = "Pthreads";
	} else {
		$pthreadsName = "NoPthreads";
	}
	return $pthreadsName;
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
	my $modelName = "UNKNOWN";

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
	my $operatorsName = "UNKNOWN";

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

