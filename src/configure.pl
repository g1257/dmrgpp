#!/usr/bin/perl 
=pod
// BEGIN LICENSE BLOCK
Copyright © 2009 , UT-Battelle, LLC
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
my ($model,$linSize,$geometry,$legOfLadder,$modelLocation,$modelFullName);
my ($geometryArgs,$nthreads);
my $wftransformation="WaveFunctionTransformation";
my $stackTemplate="MemoryStack";
my $createInput="y";
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

createInput() if ($createInput=~/^y/i);

createObserverDriver();


sub welcome
{
	print "This is DMRG++ $brand\n";
	print "This script will help you create a Makefile and main driver for your program\n";
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
	print "Available: Hubbard, Heisenberg or FeBasedSc or TjOneOrbital\n";
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
	
	$modelLocation="Models/HubbardOneBand";
	$modelLocation="Models/HeisenbergSpinOneHalf" if ($model=~/Heisenberg/i); 
	$modelLocation="Models/FeAsModel" if ($model=~/febasedsc/i);
	$modelLocation="Models/TjOneOrbital" if ($model=~/tjoneorbital/i);
	
	if ($model=~/hubbard/i or $model=~/febasedsc/i or $model=~/tjoneorbital/i) {
		$connectors="hoppings";
	} elsif ($model=~/heisenberg/i) {
		$connectors="jvalues";
	}
	
	$connectors2="jvalues" if ($model=~/tjoneorbital/i);

	if ($model=~/febasedsc/i) {
		$geometry="ladderfeas";
	} else {
		print "What geometry do you want to use?\n";
		print "Available: 1D, Ladder or LadderFeAs\n";
		print "Default is: 1D (press ENTER): ";
		$_=<STDIN>;
		chomp;
		if ($_ eq "" or $_ eq "\n") {
			$_="1d";
		}
		$geometry=$_;
	}
	
	print "Please enter the linker flags, LDFLAGS\n";
	print "Available: Any\n";
	print "Default is: $lapack (press ENTER): ";
	$_=<STDIN>;
	chomp;
	if ($_ eq "" or $_ eq "\n") {
		$_=" $lapack ";
	}
	$lapack = $_;
	
	if ($geometry=~/ladder$/i) {
		print "Enter the leg of ladder\n";
		print "Available: any even number\n";
		print "Default is: 2 (press ENTER): ";
		$_=<STDIN>;
		chomp;
		if ($_ eq "" or $_ eq "\n") {
			$_=2;
		}
		$legOfLadder=$_;
	}
	
	print "What targetting do you want?\n";
	print "Available: GroundStateTargetting TimeStepTargetting\n";
	print "Default is: GroundStateTargetting (press ENTER): ";
	$_=<STDIN>;
	chomp;
	if ($_ eq "" or $_ eq "\n") {
		$_="GroundStateTargetting";
	}
	$targetting = $_;
	
	print "Do you want to create an input file?\n";
	print "Available: y or n\n";
	print "Default is: y (press ENTER): ";
	$_=<STDIN>;
	chomp;
	if ($_ eq "" or $_ eq "\n") {
		$_="y";
	}
	$createInput = $_;
	return unless ($_=~/^y/i);
	
	print "Enter the number of kept states for the infinite loop\n";
	print "Available: any\n";
	print "Default is: 64 (press ENTER): ";
	$_=<STDIN>;
	chomp;
	if ($_ eq "" or $_ eq "\n") {
		$_=64;
	}
	$infiniteKeptStates=$_;

	print "Enter the total number of sites (system+environment)\n";
	print "Available: any\n";
	print "Default is: 16 (press ENTER): ";
	$_=<STDIN>;
	chomp;
	if ($_ eq "" or $_ eq "\n") {
		$_=16;
	}
	$linSize=$_;

	askAboutFiniteLoops();		
	
	print "Enter the value of the $connectors\n";
	print "Available: any\n";
	print "Default is: 1.0 (press ENTER): ";
	$_=<STDIN>;
	chomp;
	if ($_ eq "" or $_ eq "\n") {
		$_="1.0";
	}
	$connectorValue=$_;
	
	if (defined($connectors2)) {
		print "Enter the value of the $connectors2\n";
		print "Available: any\n";
		print "Default is: 1.0 (press ENTER): ";
		$_=<STDIN>;
		chomp;
		if ($_ eq "" or $_ eq "\n") {
			$_="1.0";
		}
		$connectorValue2=$_;
	}
	
	if ($model=~/hubbard/i) {
		askQuestionsHubbard();
	} elsif ($model=~/febasedsc/i) {
		askQuestionsFeAs();
	} elsif ($model=~/tjoneorbital/i) {
#askQuestionsTjOneOrbital(); # FIXME: ask questions about potential
	}

	print "Do you want to run with SU(2) symmetry enabled?\n";
	print "Available: yes or no\n";
	print "Default is: no (press ENTER): ";
	$_=<STDIN>;
	chomp;
	if ($_ eq "" or $_ eq "\n") {
		$_="no";
	}
	$su2Symmetry=$_;
	if ($su2Symmetry=~/y/i) {
		print "There will be ".($linSize)." sites in total, so...\n";
		print "how many total electrons should I consider?\n";
		print "Default is: ".($linSize)." (press ENTER): ";
		$_=<STDIN>;
		chomp;
		if ($_ eq "" or $_ eq "\n") {
                	$_=$linSize;
        	}
		$electrons=$_;
		print "If there are ".($linSize)." sites in total...\n";
                print "What value of angular momentum j do you want?\n";
                print "Default is: 0 (press ENTER): ";
                $_=<STDIN>;
                chomp;
                if ($_ eq "" or $_ eq "\n") {
                        $_=0;
                }
                $momentumJ=$_;
	}


}
	
sub askQuestionsHubbard()
{ 	
	print "Enter the value of the Hubbard U values\n";
	print "Available: any\n";
	print "Default is: 1.0 (press ENTER): ";
	$_=<STDIN>;
	chomp;
	if ($_ eq "" or $_ eq "\n") {
		$_="1.0";
	}
	$hubbardUvalue=$_;
	
	print "Enter the value of the potential values\n";
	print "Available: any\n";
	print "Default is: 0.0 (press ENTER): ";
	$_=<STDIN>;
	chomp;
	if ($_ eq "" or $_ eq "\n") {
		$_="0.0";
	}
	$potentialVvalue=$_;
}

sub askQuestionsFeAs()
{
	print STDERR "Please check your input file is written correctly\n";
}

sub createMakefile
{
	system("cp Makefile Makefile.bak") if (-r "Makefile");
	my $compiler = compilerName();
	my $headerFiles = join(' ', glob("Engine/*.h Models/*/*.h Geometries/*.h"));
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


sub createInput
{
	system("cp input.inp input.bak") if (-r "input.inp");
	open(FOUT,">input.inp") or die "Cannot open file input.inp for writing: $!\n";
	my $connectorValues=createConnectors($connectorValue);
	my $version=getVersion();
	my $qns = "2 0.5 0.5";
	my $inputForTimeEvolution = getTimeEvolutionInput();
	my $geometryName = getGeometryName();
	
	if ($su2Symmetry=~/y/i) {
		$_ = $electrons/(2*$linSize);
		$momentumJ /= ($linSize);
		$qns="3 $_ $_ $momentumJ\n";
		$su2Symmetry="useSu2Symmetry";
	} else {
		$su2Symmetry="nosu2";
	}
	my $terms = 1; # it is 2 for t-j model
	my $edof = 1;
	print FOUT "TotalNumberOfSites=$linSize\n";
	print FOUT "NumberOfTerms=$terms\n";
	print FOUT "DegreesOfFreedom=$edof\n";
	print FOUT "GeometryKind=$geometryName\n";
	print FOUT "GeometryOptions=ConstantValues\n";
	print FOUT "Connectors 1 1.0\n"; # FIXME only valid for hubbard model
	if ($model=~/febasedsc/i) {
	$qns = "3 1.0 1.0 0.0"; 
	print FOUT<<EOF;
hoppings 2 2
-0.058 0
0 -0.2196
hoppings 2 2
-0.2196 0
0 -0.058
hoppings 2 2
+0.20828 +0.079
+0.079 +0.20828
hoppings 2 2
+0.20828 -0.079
-0.079 +0.20828
EOF
	} else {
		# FIXME
	}
	
	if ($model=~/hubbard/i or $model=~/febasedsc/i) {	
		my $hubbardu=createHubbardU();
		print FOUT<<EOF;
hubbardU	$hubbardu
EOF
	}
	if ($model=~/hubbard/i or $model=~/febasedsc/i or $model=~/tjoneorbital/i) {
		my $potentialv=createPotentialV();

print FOUT<<EOF;
potentialV	$potentialv
density=1.0
EOF
	}

	my $hasThreads="";
	$hasThreads = "hasThreads" if ($pthreads);

	print FOUT<<EOF;
SolverOptions=hasQuantumNumbers,nowft,$su2Symmetry,$hasLoops,$hasThreads,$targetting
Version=$version
OutputFile=data.txt
InfiniteLoopKeptStates=$infiniteKeptStates
FiniteLoops $finiteLoops
TargetQuantumNumbers $qns
   
EOF
	print FOUT "Threads=$nthreads\n" if ($pthreads);
	print FOUT "$inputForTimeEvolution\n\n" if ($targetting=~/timestep/i);
	print STDERR "File input.inp has been written\n";
	close(FOUT);
}

sub getTimeEvolutionInput
{
	return "#NO_TIME_EVOLUTION" unless ($targetting=~/timestep/i);
	my $ret = <<EOF;
FILENAME tst.txt
TIMESTEP 0.1
MAXTIMES 4 
ADVANCEEACH 5 
SITE  2 10 11 
STARTINGLOOP 2 0 0 

TIMEEVOLUTION raw 
RAW_MATRIX 4 4
0.0    0.0    0.0   0.0
0.0   0.0    0.0   -1.0
0.0    0.0    0.0   1.0 
0.0    0.0    0.0   0.0 
FERMIONSIGN -1
JMVALUES 0 0
angularFactor 1

TIMEEVOLUTION raw
RAW_MATRIX 4 4
0.0    0.0    0.0   0.0
1.0    0.0    0.0   0.0
1.0    0.0    0.0   0.0
0.0    0.0    0.0   0.0
FERMIONSIGN -1
JMVALUES 0 0
angularFactor 1
EOF
	return $ret;
}

sub createConnectors
{
	my ($value)=@_;
	if ($geometry=~/ladder$/i) {
		return createConnectorsLadder($value);
	} elsif ($geometry=~/1d/i) {
		return createConnectorsChain($value);
	}
	return "UNDEFINED_GEOMETRY"; # should never reach here
}

sub createConnectorsChain
{
	my ($value)=@_;
	my $sbOpen=""; #"["
	my $sbClose=""; #"]"
	my $spacer=" "; #",";
	my $ret=$sbOpen;
	my ($i,$j);
	for ($j=0;$j<$linSize;$j++) {
		$ret=$ret.$sbOpen;
		for ($i=0;$i<$linSize;$i++) {
			if (abs(int($i-$j))==1) {
				$ret=$ret.$value;
			} else {
				$ret=$ret."0.0";
			}
			$ret=$ret.$spacer unless ($i==$linSize-1);
		}
		$ret=$ret."$sbClose\n";
		$ret=$ret.$spacer if ($j<$linSize-1);
	}
	$ret=$ret.$sbClose unless ($sbClose eq "");
		
	return $ret;
}

sub createConnectorsLadder
{
	my ($value)=@_;
	my ($x,$y,$i,$j);
	my @matrix;
	for ($x=0;$x<$linSize;$x++) {
		for ($y=0;$y<$linSize;$y++) {
			$matrix[$x][$y]=0;
		}
	}
	
	for ($x=0;$x<int($linSize/$legOfLadder);$x++) {
		for ($y=0;$y<$legOfLadder;$y++) {
			$i = $y +$x*$legOfLadder;
			if ($y+1<$legOfLadder) {
				$j=$y+1 + $x*$legOfLadder;
				$matrix[$i][$j]=$matrix[$j][$i]=$value;
			}
			if ($x+1<int($linSize/$legOfLadder)) {
				$j=$y + ($x+1)*$legOfLadder;
				$matrix[$i][$j]=$matrix[$j][$i]=$value;
			}
		}
	}
	my $sbOpen=""; #"["
	my $sbClose=""; #"]"
	my $spacer=" "; #",";
	my $ret=$sbOpen;
	for ($x=0;$x<$linSize;$x++) {
		$ret=$ret.$sbOpen;
		for ($y=0;$y<$linSize;$y++) {
			$ret=$ret.$matrix[$x][$y];
			$ret=$ret.$spacer unless ($y==$linSize-1);
		}
		$ret=$ret."$sbClose\n";
		$ret=$ret.$spacer unless ($x==$linSize-1);
	}
	$ret=$ret.$sbClose unless ($sbClose eq "");
	
	return $ret;
}	
	
sub createHubbardU
{
	#my $ret="[";
	return "4 0.0 0.0 0.0 0.0" if ($model=~/febasedsc/i);
	my $ret="$linSize ";
	my ($j);
	for ($j=0;$j<$linSize-1;$j++) {
		$ret=$ret."$hubbardUvalue ";
	}
	$ret=$ret."$hubbardUvalue"; #]";
	return $ret;
}

sub createPotentialV
{
	#my $ret="[";
	return "4 0.0 0.0 0.0 0.0" if ($model=~/febasedsc/i);
	my $ret=" ".(2*$linSize)." ";
	my ($j);
	for ($j=0;$j<2*$linSize-1;$j++) {
		$ret=$ret."$potentialVvalue ";
	}
	$ret=$ret."$potentialVvalue"; #]";
	return $ret;
}

sub getVersion
{
	my $version="NOGIT";
	my $hasGit=0;
	# First try to see if git is available
	my $tmp=system("git log -1 >& /dev/null");
	if ($tmp==0) {
		$hasGit=1;
	}
	$version=getGitVersion() if ($hasGit);
	
	 
	#If not then read it from a file
	if ($version=~/NOGIT/i) {
		if (open(VERSIONF,"version.txt")) {
			while(<VERSIONF>) {
				chomp;
				if (/commit (.*$)/) {
					$version=$1;
				}
			}
			close(VERSIONF);
		} else {
			$version="UNDEFINED";
			print STDERR "$0: (WARNING): Could not open version.txt for reading: $!\n";
		}
	}
	return $version;
}


sub getGitVersion
{
	my $version="NOGIT";
	open(PIPE,"git log -1 |") or return $version;
	while(<PIPE>) {
		chomp;
		if (/commit (.*$)/) {
			$version=$1;
		}
	}
	close(PIPE);
	if (open(VERSIONF,">version.txt")) {
		print VERSIONF "commit $version\n";
		close(VERSIONF);
	} else {
		print STDERR "$0: (WARNING): Could not open version.txt for writing: $!\n";
	}
	return $version;

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
		SparseMatrixType matrixN(model.getOperator("c",0));
		SparseMatrixType matrix2;
		transposeConjugate(matrix2,matrixN);
		SparseMatrixType A;
		multiply(A,matrix2,matrixN);
		Su2RelatedType su2Related1;
		OperatorType opN(A,1,std::pair<size_t,size_t>(0,0),1,su2Related1);
		std::cout<<"site value time\\n";
		for (size_t i0 = 0;i0<observe.size();i0++) {
			FieldType tmp = observe.template onePoint<ApplyOperatorType>(i0,opN);
			std::cout<<observe.site()<<" "<<tmp<<" "<<observe.time()<<"\\n";
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
	if  ($modelName=~/heisenberg/i) {
		print print OBSOUT<<EOF;
		// Si^+ Sj^-
		const psimag::Matrix<FieldType>& sPlus = model.getOperator("+");
		const psimag::Matrix<FieldType>& v4=observe.correlations(n,sPlus,sPlus,1);
		if (concurrency.root()) {
			std::cout<<"OperatorSplus:\\n";
			std::cout<<v4;
		}
	
		// Si^- Sj^+
		const psimag::Matrix<FieldType>& sMinus = model.getOperator("-");
		const psimag::Matrix<FieldType>& v5=observe.correlations(n,sMinus,sMinus,1);
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

sub askAboutFiniteLoops
{
	print "Do you want to do finite loops?\n";
        print "Available: y or n\n";
        print "Default is: y (press ENTER): ";
        $_=<STDIN>;
        chomp;
        if ($_ eq "" or $_ eq "\n") {
                $_="y";
        }
	$finiteLoops = "";
	$hasLoops="";
	if ($_=~/n/i) {
		$finiteLoops="1 1 100 0";
		$hasLoops="nofiniteloops";
		return;
	}
	my ($log,$m);
	my $x = $linSize/2-1;
	my $counter=1;
	while (1) {
		$m = addFiniteLoop();
		print "Do you want to add another finite loop?\n";
        	print "Available: y or n\n";
        	print "Default is: n (press ENTER): ";
        	$_=<STDIN>;
        	chomp;
        	if ($_ eq "" or $_ eq "\n") {
                	$_="n";
        	}
		$log = 0;
		last if ($_=~/n/i);
		$log =0;
		$finiteLoops = $finiteLoops." $x $m $log -$x $m $log -$x $m $log $x $m $log ";
		$counter++;
	}
	$finiteLoops = (4*$counter)." ".$finiteLoops." $x $m $log -$x $m $log -$x $m $log $x $m $log ";
}

sub addFiniteLoop
{
	print "What's the m for this loop?\n";
	print "Available: any\n";
	print "Default is: 60 (press ENTER): ";
	$_=<STDIN>;
	chomp;
	if ($_ eq "" or $_ eq "\n") {
		$_=60;
	}
	return $_;

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

sub getGeometryName
{
	my $geometryName = "UNKNOWN";
	if ($geometry=~/1d/i) {
		$geometryName = "chain";
	} elsif ($geometry=~/ladder$/i) {
		$geometryName = "ladder";
	} elsif ($geometry=~/ladderfeas/i) {
		$geometryName = "ladderx";
	}
	return $geometryName;
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

	if ($model=~/hubbard/i) {
		$modelName = "ModelHubbard";
	} elsif ($model=~/heisenberg/i) {
		$modelName="ModelHeisenberg";
	} elsif ($model=~/febasedsc/i) {
		$modelName = "ModelFeBasedSc";
	} elsif ($model=~/tjoneorbital/i) {
		$modelName = "TjOneOrbital";
	}
	return $modelName;
}

sub getOperatorsName()
{
	my $operatorsName = "UNKNOWN";

	if ($model=~/hubbard/i) {
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
	my @tryThis = ("g++4","g++");
	my $ret;
	my $compiler;
	foreach my $comp (@tryThis) {
		system("$comp >& /dev/null");
		next if ($? == -1 or ($? & 127) );
		$ret = $? >> 8;
		if ($ret < 2) {
			$compiler = $comp;
			last;
		} 
			
	}
	return $compiler if defined $compiler; 
	die "$0: No suitable compiler found\n";
}

