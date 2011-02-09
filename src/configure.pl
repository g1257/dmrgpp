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
$DynamicTargetting = "DynamicTargettingEmpty" if ($hasGsl=~/n/i);

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
	print "Available: Hubbard, Heisenberg or FeBasedSc or TjOneOrbital or ExtendedHubbard1Orbital\n";
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
	$modelFullName="FeAsModel" if ($model=~/tjoneorbital/i);
	$modelFullName="ExtendedHubbard1Orb" if ($model=~/extendedhubbard1orb/i);
	
	$modelLocation="Models/HubbardOneBand";
	$modelLocation="Models/HeisenbergSpinOneHalf" if ($model=~/Heisenberg/i); 
	$modelLocation="Models/FeAsModel" if ($model=~/febasedsc/i);
	$modelLocation="Models/TjOneOrbital" if ($model=~/tjoneorbital/i);
	$modelLocation="Models/ExtendedHubbard1Orb" if ($model=~/extendedhubbard1orb/i);
	
	if ($model=~/hubbard/i or $model=~/febasedsc/i or $model=~/tjoneorbital/i
		or $model=~/extendedhubbard1orb/i) {
		$connectors="hoppings";
	} elsif ($model=~/heisenberg/i) {
		$connectors="jvalues";
	}
	
	$connectors2="jvalues" if ($model=~/tjoneorbital/i);
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
	if ($modelLocation=~/extendedhubbard1orb/i) {
		$modelLocation = $modelLocation." -IModels/HubbardOneBand ";
	}
	my $litProgTargets = getLitProgTargets(\@litProgFiles);
	open(FOUT,">Makefile") or die "Cannot open Makefile for writing: $!\n";
print FOUT<<EOF;
# DO NOT EDIT!!! Changes will be lost. Modify configure.pl instead
# This Makefile was written by configure.pl
# DMRG++ ($brand) by G.A.
# Platform: $platform
# MPI: $mpi

LDFLAGS =    $lapack  $gslLibs $pthreadsLib
CPPFLAGS = -Werror -Wall  -IEngine -I$modelLocation -IGeometries -I$PsimagLite
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
#include "DensityMatrix.h"
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

typedef double MatrixElementType;
typedef std::complex<MatrixElementType> ComplexType;
typedef  Dmrg::CrsMatrix<ComplexType> MySparseMatrixComplex;
typedef  Dmrg::CrsMatrix<MatrixElementType> MySparseMatrixReal;

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
	typedef ModelHelperTemplate<OperatorsType,ReflectionSymmetryType,MyConcurrency> ModelHelperType;
	typedef ModelTemplate<ModelHelperType,MySparseMatrix,GeometryType,PsimagLite::$pthreadsName> ModelType;
	
	typedef DmrgSolver<
			InternalProductTemplate,
			DensityMatrix,
			ModelType,
			MyConcurrency,
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
#include "DensityMatrix.h" // only used for types
#include "TimeStepTargetting.h" // only used for types
#include "GroundStateTargetting.h" // only used for types

using namespace Dmrg;

typedef double RealType;
typedef std::complex<RealType> ComplexType;
$chooseRealOrComplexForObservables

typedef PsimagLite::$concurrencyName<RealType> MyConcurrency;
typedef PsimagLite::IoSimple::In IoInputType;

template<typename ModelType,typename ObserverType,typename SparseMatrixType,
typename OperatorType,typename TargettingType>
void measureTimeObs(const ModelType& model,ObserverType& observe, size_t numberOfSites)
{
	typedef typename OperatorType::Su2RelatedType Su2RelatedType;
	typedef typename TargettingType::ApplyOperatorType ApplyOperatorType;
	SparseMatrixType matrixNup(model.getOperator("nup"));
	SparseMatrixType matrixNdown(model.getOperator("ndown"));
	SparseMatrixType A;
	Su2RelatedType su2Related1;
	A.makeDiagonal(4,1.0);
	OperatorType opIdentity(A,1,std::pair<size_t,size_t>(0,0),1,su2Related1);
	FieldType superDensity = observe.template onePoint<ApplyOperatorType>(0,opIdentity);
	std::cout<<"SuperDensity(Weight of the timeVector)="<<superDensity<<"\\n";

	multiply(A,matrixNup,matrixNdown);
	std::cout<<"#Using Matrix A:\\n";
	for (size_t i=0;i<A.rank();i++) {
		for (size_t j=0;j<A.rank();j++)
			std::cout<<"#A("<<i<<","<<j<<")="<<A(i,j)<<" ";
		std::cout<<"\\n";
	}
	OperatorType opN(A,1,std::pair<size_t,size_t>(0,0),1,su2Related1);
	std::cout<<"site nupNdown(gs) nupNdown(timevector) time\\n";
	for (size_t i0 = 0;i0<observe.size();i0++) {
		// for g.s. use this one:
		observe.setBrackets("gs","gs");
		FieldType tmp1 = observe.template onePoint<ApplyOperatorType>(i0,opN);
		// for time vector use this one:
		observe.setBrackets("time","time");
		FieldType tmp2 = observe.template onePoint<ApplyOperatorType>(i0,opN);
		std::cout<<observe.site()<<" "<<tmp1<<" "<<tmp2<<" "<<observe.time()<<"\\n";
		if (observe.isAtCorner()) { // also calculate next or prev. site:
			size_t x = (observe.site()==1) ? 0 : numberOfSites-1;
			// do the corner case
			// for g.s. use this one:
			observe.setBrackets("gs","gs");
			FieldType tmp1 = observe.template onePoint<ApplyOperatorType>(i0,opN,true);
			// for time vector use this one:
			observe.setBrackets("time","time");
			FieldType tmp2 = observe.template onePoint<ApplyOperatorType>(i0,opN,true);
			std::cout<<x<<" "<<tmp1<<" "<<tmp2<<" "<<observe.time()<<"\\n";
		}
	}
	// measuring charge:
	A = matrixNup;
	A += matrixNdown;
	std::cout<<"#Using Matrix A:\\n";
	for (size_t i=0;i<A.rank();i++) {
		for (size_t j=0;j<A.rank();j++)
			std::cout<<"#A("<<i<<","<<j<<")="<<A(i,j)<<" ";
		std::cout<<"\\n";
	}
	OperatorType opCharge(A,1,std::pair<size_t,size_t>(0,0),1,su2Related1);
	std::cout<<"site nUp+nDown(gs) nup+ndown(timevector) time\\n";
	for (size_t i0 = 0;i0<observe.size();i0++) {
		//for g.s. use this one:
		observe.setBrackets("gs","gs");
		FieldType tmp1 = observe.template onePoint<ApplyOperatorType>(i0,opCharge);
		// for h.d. use this one:
		observe.setBrackets("time","time");
		FieldType tmp2 = observe.template onePoint<ApplyOperatorType>(i0,opCharge);
		std::cout<<observe.site()<<" "<<tmp1<<" "<<tmp2<<" "<<observe.time()<<"\\n";
		if (observe.isAtCorner()) { // also calculate next or prev. site:
			size_t x = (observe.site()==1) ? 0 : numberOfSites-1;
			// do the corner case
			// for g.s. use this one:
			observe.setBrackets("gs","gs");
			FieldType tmp1 = observe.template onePoint<ApplyOperatorType>(i0,opCharge,true);
			// for time vector use this one:
			observe.setBrackets("time","time");
			FieldType tmp2 = observe.template onePoint<ApplyOperatorType>(i0,opCharge,true);
			std::cout<<x<<" "<<tmp1<<" "<<tmp2<<" "<<observe.time()<<"\\n";
		}
	}

}

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
	ObserverType observe(io,n-2,hasTimeEvolution,model,concurrency,verbose);
	
	if (hasTimeEvolution) {
		measureTimeObs<ModelType,ObserverType,SparseMatrixType,
			OperatorType,TargettingType>(model,observe,geometry.numberOfSites());
		return observe.endOfData(); // return here for testing only 
	}
	observe.setBrackets("gs","gs");
EOF
	if  ($modelName=~/heisenberg/i) {
	} else {
print OBSOUT<<EOF;
	// OPERATOR C:
	if (obsOptions.find("cc")!=std::string::npos) {
		PsimagLite::Matrix<RealType> opC1 = model.getOperator("c",0,0); // c_{0,0} spin up
		PsimagLite::Matrix<FieldType> opC = opC1;	
		PsimagLite::Matrix<FieldType> opCtranspose = utils::transposeConjugate(opC1);
		const PsimagLite::Matrix<FieldType>& v=observe.correlations(n,opC,opCtranspose,-1);;
		if (concurrency.root()) {
			std::cout<<"OperatorC:\\n";
			std::cout<<v;
		}
	}
	
	// OPERATOR N
	if (obsOptions.find("nn")!=std::string::npos) {
		PsimagLite::Matrix<RealType> opN1 = model.getOperator("n");
		PsimagLite::Matrix<FieldType> opN = opN1;
		const PsimagLite::Matrix<FieldType>& v2=observe.correlations(n,opN,opN,1);
		if (concurrency.root()) {
			std::cout<<"OperatorN:\\n";
			std::cout<<v2;
		}
	}
EOF
	}
print OBSOUT<<EOF;
	// OPERATOR SZ
	PsimagLite::Matrix<FieldType>* v3;
	if (obsOptions.find("szsz")!=std::string::npos) {
		PsimagLite::Matrix<RealType> Sz1 = model.getOperator("z");
		PsimagLite::Matrix<FieldType> Sz = Sz1;
		v3=new PsimagLite::Matrix<FieldType>(observe.correlations(n,Sz,Sz,1));
		if (concurrency.root()) {
			std::cout<<"OperatorSz:\\n";
			std::cout<<(*v3);
		}
	}
EOF
	if  ($modelName=~/febasedsc/i) {
print print OBSOUT<<EOF;
	// TWO-POINT DELTA-DELTA^DAGGER:
	if (obsOptions.find("dd")!=std::string::npos && 
			geometry.label(0).find("ladder")==std::string::npos) {
		const PsimagLite::Matrix<FieldType>& oDelta = model.getOperator("d");
		PsimagLite::Matrix<FieldType> oDeltaT;
		utils::transposeConjugate(oDeltaT,oDelta);
		const PsimagLite::Matrix<FieldType>& vDelta=observe.correlations(n,oDelta,oDeltaT,1);
		if (concurrency.root()) {
			std::cout<<"DeltaDeltaDagger:\\n";
			std::cout<<vDelta;
		}
	}
	

	// FOUR-POINT DELTA-DELTA^DAGGER:
	if (obsOptions.find("dd4")!=std::string::npos &&
		geometry.label(0).find("ladder")!=std::string::npos) {
		for (size_t g=0;g<16;g++) {
			std::vector<FieldType> fpd;
			std::vector<size_t> gammas(4,0); // TODO: WRITE COMBINATIONS
			gammas[0] = (g & 1);
			gammas[1] = (g & 2)>>1;
			gammas[2] = (g & 4) >> 2;
			gammas[3] = (g & 8) >> 3;
			observe.fourPointDeltas(fpd,geometry.numberOfSites(),gammas,model);
			for (size_t step=0;step<fpd.size();step++) {
				// step --> (i,j) FIXME
				std::cout<<step<<" "<<g<<" "<<fpd[step]<<"\\n";
			}
		}
	} // if dd4
EOF
	}
	if  ($modelName=~/heisenberg/i) {
	print print OBSOUT<<EOF;
		// Heisenberg model: needs formatting and also obsOptions: FIXME
		// Si^+ Sj^-
		const PsimagLite::Matrix<FieldType>& sPlus = model.getOperator("+");
		PsimagLite::Matrix<FieldType> sPlusT = utils::transposeConjugate(sPlus);
		const PsimagLite::Matrix<FieldType>& v4=observe.correlations(n,sPlus,sPlusT,1);
		if (concurrency.root()) {
			std::cout<<"OperatorSplus:\\n";
			std::cout<<v4;
		}
	
		// Si^- Sj^+
		const PsimagLite::Matrix<FieldType>& sMinus = model.getOperator("-");
		PsimagLite::Matrix<FieldType> sMinusT = utils::transposeConjugate(sMinus);
		const PsimagLite::Matrix<FieldType>& v5=observe.correlations(n,sMinus,sMinusT,1);
		if (concurrency.root()) {
			std::cout<<"OperatorMinus:\\n";
			std::cout<<v5;
		}
	
		PsimagLite::Matrix<FieldType> spinTotal(v5.n_row(),v5.n_col());
		
		for (size_t i=0;i<spinTotal.n_row();i++) 
			for (size_t j=0;j<spinTotal.n_col();j++) spinTotal(i,j) = 0.5*(v4(i,j) +v5(i,j)) + (*v3)(i,j);
	
		if (concurrency.root()) {
			std::cout<<"SpinTotal:\\n";
			std::cout<<spinTotal;
		}
EOF
	}
	print OBSOUT<<EOF;
	return observe.endOfData();
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
	typedef ModelHelperTemplate<OperatorsType,ReflectionSymmetryType,ConcurrencyType> ModelHelperType;
	typedef ModelTemplate<ModelHelperType,MySparseMatrix,GeometryType,PsimagLite::$pthreadsName> ModelType;
	
	typedef DmrgSolver<
                        InternalProductTemplate,
                        DensityMatrix,
                        ModelType,
                        MyConcurrency,
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
	std::string options = "cc,nn,szsz";
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
	typedef CrsMatrix<RealType> MySparseMatrixReal;
	if (hasTimeEvolution)
		mainLoop<ParametersModelType,GeometryType,MyConcurrency,IoInputType,$modelName,
			ModelHelperLocal,InternalProductOnTheFly,VectorWithOffsets,
			$targetting,MySparseMatrixReal>(mp,geometry,hasTimeEvolution,concurrency,io,dmrgSolverParams.filename,options);
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

