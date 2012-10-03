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
my ($pthreads,$pthreadsLib)=(0,"");
my $brand= "v2.0";

my $gslLibs = " -lgsl  -lgslcblas ";
$gslLibs =" " if ($hasGsl=~/n/i);

system("make clean");

guessPlatform();

welcome();

askQuestions();

createMakefile();

#createDriver();

#createObserverDriver();


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
	unlink("Engine/Version.h");
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
CPPFLAGS = -Werror -Wall  -IEngine -IModels/HubbardOneBand -IModels/HeisenbergSpinOneHalf -IModels/ExtendedHubbard1Orb  -IModels/FeAsModel -IModels/FeAsBasedScExtended -IModels/Immm  -I$PsimagLite -I$PsimagLite/Geometry $usePthreadsOrNot
EOF
if ($mpi) {
	print FOUT "CXX = mpicxx -O3 -DNDEBUG \n";
} else {
	print FOUT "CXX = $compiler  -O3 -DNDEBUG\n";
	print FOUT "#Comment out line below for debugging (COMMENT ALSO THE strip commands): \n";
	print FOUT "#CXX = $compiler -g3 #ALSO COMMENT OUT strip command below\n";
}
print FOUT<<EOF;
EXENAME = dmrg
all: \$(EXENAME)

dmrg:  dmrg.o gitrev
	\$(CXX) -o dmrg dmrg.o \$(LDFLAGS)  
	strip dmrg

correctionVectorMulti: correctionVectorMulti.o
	\$(CXX) -o correctionVectorMulti correctionVectorMulti.o \$(LDFLAGS)

observe:  observe.o
	\$(CXX) -o observe observe.o \$(LDFLAGS)
	strip observe

# dependencies brought about by Makefile.dep
%.o: %.cpp Makefile gitrev Engine/Version.h
	\$(CXX) \$(CPPFLAGS) -c \$< 

Makefile.dep: Engine/Version.h dmrg.cpp
	\$(CXX) \$(CPPFLAGS) -MM dmrg.cpp  > Makefile.dep

Engine/Version.h: gitrev
	./gitrev > Engine/Version.h	

gitrev: gitrev.o
	\$(CXX) -o gitrev gitrev.o \$(LDFLAGS)

gitrev.o:
	\$(CXX) \$(CPPFLAGS) -c gitrev.cpp 

doc:
	cd ../doc; make

# dependencies brought about by MakefileObserver.dep
observe.o: Engine/Version.h
	\$(CXX) \$(CPPFLAGS) -c observe.cpp

MakefileObserver.dep: Engine/Version.h observe.cpp
	\$(CXX) \$(CPPFLAGS) -MM observe.cpp  > MakefileObserver.dep

clean:
	rm -f core* \$(EXENAME) *.o *.dep Engine/Version.h gitrev

include Makefile.dep
include MakefileObserver.dep
######## End of Makefile ########

EOF
	close(FOUT);
	print STDERR "File Makefile has been written\n";
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


sub compilerName
{
	return "g++";
	my @tryThis = ("g++","g++4");
	my $ret;
	my $compiler;
	foreach my $comp (@tryThis) {
		my $ret = system("$comp > /dev/null 2>/dev/null");
		if ($ret==0) {
			$compiler = $comp;
			last;
		} 
			
	}
	return $compiler if defined $compiler; 
	die "$0: No suitable compiler found\n";
}

