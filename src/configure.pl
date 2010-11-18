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

my $mpi=0;
my $platform="linux";
my $lapack="-llapack";
my $connectors;
my $connectorValue=1.0;
my $hubbardUvalue=1.0;
my $potentialVvalue=0.0;
my ($infiniteKeptStates,$finiteLoops,$hasLoops);
my ($model,$linSize,$modelLocation,$modelFullName);
my ($geometryArgs,$nthreads);
my $wftransformation="WaveFunctionTransformation";
my $stackTemplate="MemoryStack";
my ($electrons,$momentumJ,$su2Symmetry);
my ($pthreads,$pthreadsLib)=(0,"");
my $brand= "v2.0";
my ($connectorsArgs,$connectorsArgs2,$dof,$connectors2,$connectorValue2);
my $targetting;

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
	
	if ($pthreads) {
		print "How many threads?\n";
		print "Available: Any\n";
		print "Default is: 2 (press ENTER): ";
		$_=<STDIN>;
		chomp;
		if ($_ eq "" or $_ eq "\n") {
			$_=2;
		}
		$nthreads = $_;
	}
	
	if ($mpi!=1) {	
		print "What stack class do you want to use?\n";
		print "Available: MemoryStack or DiskStack\n";
		print "Default is: MemoryStack (press ENTER): ";
		$_=<STDIN>;
		chomp;
		if ($_ eq "" or $_ eq "\n") {
			$_="MemoryStack";
		}
		$stackTemplate=$_;
	} else {
		$stackTemplate="MemoryStack";
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
	print "Available: GroundStateTargetting TimeStepTargetting\n";
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
	if ($modelLocation=~/extendedhubbard1orb/i) {
		$modelLocation = $modelLocation." -IModels/HubbardOneBand ";
	}
	open(FOUT,">Makefile") or die "Cannot open Makefile for writing: $!\n";
print FOUT<<EOF;
# DO NOT EDIT!!! Changes will be lost. Modify configure.pl instead
# This Makefile was written by configure.pl
# DMRG++ ($brand) by G.A.
# Platform: $platform
# MPI: $mpi

LDFLAGS =    $lapack  -lm $pthreadsLib
CPPFLAGS = -Werror -Wall -I../PartialPsimag -IEngine -I$modelLocation -IGeometries
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

clean:
	rm -f core* \$(EXENAME) *.o *.ii *.tt

######## End of Makefile ########

EOF
	close(FOUT);
	print STDERR "File Makefile has been written\n";
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
#include "$wftransformation.h"
#include "DensityMatrix.h"
#include "$stackTemplate.h"
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
//typedef BlockMatrix<MatrixElementType,psimag::Matrix<MatrixElementType> > MyDenseBlockMatrix;
typedef Dmrg::$concurrencyName<MatrixElementType> MyConcurrency;
typedef $parametersName<MatrixElementType> ParametersModelType;
typedef Geometry<MatrixElementType> GeometryType;;
typedef  IoSimple MyIo;

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
	template<template<typename,typename> class,template<typename,typename> class,
		typename,typename,typename,typename,template<typename> class> class TargettingTemplate,
	typename MySparseMatrix
>
void mainLoop(ParametersModelType& mp,GeometryType& geometry,ParametersSolverType& dmrgSolverParams,
		ConcurrencyType& concurrency, IoInputType& io,bool hasTimeEvolution)
{
	typedef ReflectionSymmetryEmpty<MatrixElementType,MySparseMatrix> ReflectionSymmetryType;
	typedef Operator<MatrixElementType,MySparseMatrix> OperatorType;
	typedef Basis<MatrixElementType,MySparseMatrix> BasisType;
	typedef $operatorsName<OperatorType,BasisType> OperatorsType;
	typedef ModelHelperTemplate<OperatorsType,ReflectionSymmetryType,MyConcurrency> ModelHelperType;
	typedef ModelTemplate<ModelHelperType,MySparseMatrix,GeometryType,$pthreadsName> ModelType;
	
	typedef DmrgSolver<
			InternalProductTemplate,
			DensityMatrix,
			ModelType,
			MyConcurrency,
			MyIo,
			$wftransformation,
			$stackTemplate,
			TargettingTemplate,
			VectorWithOffsetTemplate
		> SolverType;
	
	//! Setup the Model
	ModelType model(mp,geometry);

	//! Read TimeEvolution if applicable:
	typedef typename SolverType::TargettingType TargettingType;
        typedef typename TargettingType::TargettingStructureType TargettingStructureType;
        TargettingStructureType tsp(io,model,hasTimeEvolution);

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
	typedef IoSimple::In IoInputType;
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
	bool hasTimeEvolution=false;
	if (dmrgSolverParams.options.find("TimeStepTargetting")!=std::string::npos) hasTimeEvolution=true;
	if (hasTimeEvolution && su2) throw std::runtime_error("Time Evolution and SU(2)"
 		" not supported at the same time yet (sorry!)\\n");
	if (!su2) {
		if (hasTimeEvolution) { 
			mainLoop<ParametersModelType,GeometryType,ParametersDmrgSolver<MatrixElementType>,MyConcurrency,
				IoInputType,
				$modelName,ModelHelperLocal,InternalProductOnTheFly,VectorWithOffsets,TimeStepTargetting,
				MySparseMatrixComplex>
			(mp,geometry,dmrgSolverParams,concurrency,io,hasTimeEvolution);
		} else {
			mainLoop<ParametersModelType,GeometryType,ParametersDmrgSolver<MatrixElementType>,MyConcurrency,
				IoInputType,
				$modelName,ModelHelperLocal,InternalProductOnTheFly,VectorWithOffset,GroundStateTargetting,
				MySparseMatrixReal>
			(mp,geometry,dmrgSolverParams,concurrency,io,hasTimeEvolution);
		}
	} else {
		 if (dmrgSolverParams.targetQuantumNumbers[2]>0) { 
			mainLoop<ParametersModelType,GeometryType,ParametersDmrgSolver<MatrixElementType>,MyConcurrency,
				IoInputType,
				$modelName,ModelHelperSu2,InternalProductOnTheFly,VectorWithOffsets,GroundStateTargetting,
				MySparseMatrixReal>
			(mp,geometry,dmrgSolverParams,concurrency,io,hasTimeEvolution);
		} else {
			mainLoop<ParametersModelType,GeometryType,ParametersDmrgSolver<MatrixElementType>,MyConcurrency,
				IoInputType,
				$modelName,ModelHelperSu2,InternalProductOnTheFly,VectorWithOffset,GroundStateTargetting,
				MySparseMatrixReal>
			(mp,geometry,dmrgSolverParams,concurrency,io,hasTimeEvolution);
		}
	}
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
	my $obsArg = "datafile,n,opInfo.n_row(),concurrency,verbose";
	if ($targetting=~/timestep/i) {
		$chooseRealOrComplexForObservables = "typedef ComplexType FieldType;\n";
		$obsArg = "datafile,tsp.filename,n,opInfo.n_row(),concurrency,verbose";
	}
	
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
#include "WaveFunctionTransformation.h"// only used for types
#include "$stackTemplate.h" // only used for types
#include "TimeStepTargetting.h" // only used for types
#include "GroundStateTargetting.h" // only used for types

using namespace Dmrg;

typedef double RealType;
typedef std::complex<RealType> ComplexType;
$chooseRealOrComplexForObservables

typedef Dmrg::$concurrencyName<RealType> MyConcurrency;
typedef  IoSimple MyIo;

	template<
	typename ParametersModelType,
	typename GeometryType,
	typename ConcurrencyType,
	typename IoInputType,
	template<typename,typename,typename,template<typename> class> class ModelTemplate,
	template<typename,typename,typename> class ModelHelperTemplate,
	template<typename,typename> class InternalProductTemplate,
	template<typename> class VectorWithOffsetTemplate,
	template<template<typename,typename> class,template<typename,typename> class,
		typename,typename,typename,typename,template<typename> class> class TargettingTemplate,
	typename MySparseMatrix
>
void mainLoop(ParametersModelType& mp,GeometryType& geometry,bool hasTimeEvolution,
		ConcurrencyType& concurrency, IoInputType& io,const std::string& datafile)
{
	typedef ReflectionSymmetryEmpty<RealType,MySparseMatrix> ReflectionSymmetryType;
	typedef Operator<RealType,MySparseMatrix> OperatorType;
	typedef Basis<RealType,MySparseMatrix> BasisType;
	typedef $operatorsName<OperatorType,BasisType> OperatorsType;
	typedef typename OperatorType::Su2RelatedType Su2RelatedType;
	typedef typename OperatorType::SparseMatrixType SparseMatrixType;
	typedef ModelHelperTemplate<OperatorsType,ReflectionSymmetryType,ConcurrencyType> ModelHelperType;
	typedef ModelTemplate<ModelHelperType,MySparseMatrix,GeometryType,$pthreadsName> ModelType;
	
	typedef DmrgSolver<
                        InternalProductTemplate,
                        DensityMatrix,
                        ModelType,
                        MyConcurrency,
                        MyIo,
                        WaveFunctionTransformation,
                        MemoryStack,
                        TargettingTemplate,
                        VectorWithOffsetTemplate
                > SolverType; // only used for types
	
	typedef typename SolverType::TargettingType TargettingType;
	typedef typename SolverType::WaveFunctionTransformationType WaveFunctionTransformationType;
	typedef typename TargettingType::VectorWithOffsetType VectorWithOffsetType;
	typedef typename TargettingType::ApplyOperatorType ApplyOperatorType;
	ModelType model(mp,geometry);

	 //! Read TimeEvolution if applicable:
        typedef typename SolverType::TargettingType TargettingType;
        typedef typename TargettingType::TargettingStructureType TargettingStructureType;
	TargettingStructureType tsp(io,model,hasTimeEvolution);
	
	size_t n=geometry.numberOfSites()/2;
	const psimag::Matrix<FieldType>& opInfo = model.getOperator("i",0,0);
	bool verbose = false;
	Observer<FieldType,VectorWithOffsetType,BasisType,IoSimple,ConcurrencyType> observe($obsArg);
	if (hasTimeEvolution) {
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
		for (int i=0;i<A.rank();i++) {
			for (int j=0;j<A.rank();j++)
				std::cout<<"#A("<<i<<","<<j<<")="<<A(i,j)<<" ";
			std::cout<<"\\n";
		}
		OperatorType opN(A,1,std::pair<size_t,size_t>(0,0),1,su2Related1);
		std::cout<<"site nupNdown(gs) nupNdown(timevector) time\\n";
		for (size_t i0 = 0;i0<observe.size();i0++) {
			// for g.s. use this one:
			FieldType tmp1 = observe.template onePointGs<ApplyOperatorType>(i0,opN);
			// for time vector use this one:
			FieldType tmp2 = observe.template onePoint<ApplyOperatorType>(i0,opN);
			std::cout<<observe.site()<<" "<<tmp1<<" "<<tmp2<<" "<<observe.time()<<"\\n";
		}
		// measuring charge:
		A = matrixNup;
		A += matrixNdown;
		std::cout<<"#Using Matrix A:\\n";
		for (int i=0;i<A.rank();i++) {
			for (int j=0;j<A.rank();j++)
				std::cout<<"#A("<<i<<","<<j<<")="<<A(i,j)<<" ";
			std::cout<<"\\n";
		}
		OperatorType opCharge(A,1,std::pair<size_t,size_t>(0,0),1,su2Related1);
		std::cout<<"site nUp+nDown(gs) nup+ndown(timevector) time\\n";
		for (size_t i0 = 0;i0<observe.size();i0++) {
			//for g.s. use this one:
			FieldType tmp1 = observe.template onePointGs<ApplyOperatorType>(i0,opCharge);
			// for h.d. use this one:
			FieldType tmp2 = observe.template onePoint<ApplyOperatorType>(i0,opCharge);
			std::cout<<observe.site()<<" "<<tmp1<<" "<<tmp2<<" "<<observe.time()<<"\\n";
		}

		return;
	}
EOF
	if  ($modelName=~/heisenberg/i) {
	} else {
		print OBSOUT<<EOF;
		const psimag::Matrix<FieldType>& opC = model.getOperator("c",0,0);
		psimag::Matrix<FieldType> opCtranspose = transposeConjugate(opC);
		const psimag::Matrix<FieldType>& v=observe.correlations(n,opC,opCtranspose,-1);;
		if (concurrency.root()) {
			std::cout<<"OperatorC:\\n";
			std::cout<<v;
		}
	
		const psimag::Matrix<FieldType>& opN = model.getOperator("n");
		const psimag::Matrix<FieldType>& v2=observe.correlations(n,opN,opN,1);
		if (concurrency.root()) {
			std::cout<<"OperatorN:\\n";
			std::cout<<v2;
		}
EOF
	}
	print OBSOUT<<EOF;
	
	const psimag::Matrix<FieldType>& Sz = model.getOperator("z");
	const psimag::Matrix<FieldType>& v3=observe.correlations(n,Sz,Sz,1);
	if (concurrency.root()) {
		std::cout<<"OperatorSz:\\n";
		std::cout<<v3;
	}
EOF
	if  ($modelName=~/febasedsc/i) {
		print print OBSOUT<<EOF;
		const psimag::Matrix<FieldType>& oDelta = model.getOperator("d");
		psimag::Matrix<FieldType> oDeltaT;
		transposeConjugate(oDeltaT,oDelta);
		const psimag::Matrix<FieldType>& vDelta=observe.correlations(n,oDelta,oDeltaT,1);
		if (concurrency.root()) {
			std::cout<<"DeltaDeltaDagger:\\n";
			std::cout<<vDelta;
		}
EOF
	}
	if  ($modelName=~/heisenberg/i) {
		print print OBSOUT<<EOF;
		// Si^+ Sj^-
		const psimag::Matrix<FieldType>& sPlus = model.getOperator("+");
		psimag::Matrix<FieldType> sPlusT = transposeConjugate(sPlus);
		const psimag::Matrix<FieldType>& v4=observe.correlations(n,sPlus,sPlusT,1);
		if (concurrency.root()) {
			std::cout<<"OperatorSplus:\\n";
			std::cout<<v4;
		}
	
		// Si^- Sj^+
		const psimag::Matrix<FieldType>& sMinus = model.getOperator("-");
		psimag::Matrix<FieldType> sMinusT = transposeConjugate(sMinus);
		const psimag::Matrix<FieldType>& v5=observe.correlations(n,sMinus,sMinusT,1);
		if (concurrency.root()) {
			std::cout<<"OperatorMinus:\\n";
			std::cout<<v5;
		}
	
		psimag::Matrix<FieldType> spinTotal(v5.n_row(),v5.n_col());
		
		for (size_t i=0;i<spinTotal.n_row();i++) 
			for (size_t j=0;j<spinTotal.n_col();j++) spinTotal(i,j) = 0.5*(v4(i,j) +v5(i,j)) + v3(i,j);
	
		if (concurrency.root()) {
			std::cout<<"SpinTotal:\\n";
			std::cout<<spinTotal;
		}
EOF
	}
	print OBSOUT<<EOF;
}

int main(int argc,char *argv[])
{
	if (argc>2) {
		std::string s = "The observer driver takes only one argument now: \\n";
		s = s + "the name of the input file. The data file is now read from the input file.\\n";
		throw std::runtime_error(s);
	}
	using namespace Dmrg;
	
	typedef   Dmrg::$concurrencyName<RealType> MyConcurrency; 
	
	MyConcurrency concurrency(argc,argv);
	
	//Setup the Geometry
        typedef Geometry<RealType> GeometryType;
	typedef IoSimple::In IoInputType;
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
			$targetting,MySparseMatrixReal>(mp,geometry,hasTimeEvolution,concurrency,io,dmrgSolverParams.filename);
	else
		mainLoop<ParametersModelType,GeometryType,MyConcurrency,IoInputType,$modelName,
			ModelHelperLocal,InternalProductOnTheFly,VectorWithOffset,
			$targetting,MySparseMatrixReal>(mp,geometry,hasTimeEvolution,concurrency,io,dmrgSolverParams.filename);
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

