#!/usr/bin/perl
=pod
Copyright (c) 2009-2014, UT-Battelle, LLC
All rights reserved

[DMRG++, Version 3.0]

*********************************************************
THE SOFTWARE IS SUPPLIED BY THE COPYRIGHT HOLDERS AND
CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED
WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A
PARTICULAR PURPOSE ARE DISCLAIMED.

Please see full open source license included in file LICENSE.
*********************************************************

=cut
use warnings;
use strict;

use lib "../../PsimagLite/scripts";
use Make;

my $mpi=0;
my $platform="linux";
my @drivers = ("dmrg","observe");
my $lapack=Make::findLapack();
my $PsimagLite="../../PsimagLite";
my ($pthreads,$pthreadsLib)=(0,"");
my $brand= "v3.0";
my $build="production";
my $floating="";

system("make clean");

guessPlatform();

welcome();

askQuestions();

createMakefile();

sub welcome
{
	print "This is DMRG++ $brand\n";
	print "This script will help you create a Makefile for your program\n";
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
	$pthreadsLib="  ";
	if ($mpi!=1 && !($platform=~/Darwin/i)) { # cannot use both mpi and pthreads
		print "Do you want to compile with pthreads enabled?\n";
		print "Available: y or n\n";
		print "Default is: y (press ENTER): ";
		$_=<STDIN>;
		chomp;
		if ($_ eq "" or $_ eq "\n") {
			$_="y";
		}
		if ($_=~/^y/i) {
			$pthreadsLib=" -lpthread ";
			$pthreads=1;
		}
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

	print "Please enter the type of build\n";
	print "Available: production, debug, callgrind, drd\n";
	print "Default is: production (press ENTER): ";
	$_=<STDIN>;
	chomp;
	if ($_ eq "" or $_ eq "\n") {
		$_="production";
	}
	$build = $_;

	print "Please enter the type of floating point to be used\n";
	print "Available: double or float\n";
	print "Default is: double (press ENTER): ";
	$_=<STDIN>;
	chomp;
	if ($_ eq "" or $_ eq "\n") {
		$_ = "double";
	}

	if ($_ eq "float") {
		$floating = " -DUSE_FLOAT ";
	}

}

sub createMakefile
{
	unlink("Engine/Version.h");
	system("cp Makefile Makefile.bak") if (-r "Makefile");
	my $compiler = "g++";
	$compiler = " mpicxx " if ($mpi);
	my $fh;
	open($fh,">Makefile") or die "Cannot open Makefile for writing: $!\n";
	my $usePthreadsOrNot = " ";
	$usePthreadsOrNot = " -DUSE_PTHREADS " if ($pthreads);
	my $optimizations = " -O3 -DNDEBUG ";
	$optimizations = " -g3 -Wextra -D_GLIBCXX_DEBUG -D_GLIBCXX_PROFILE " if ($build eq "debug");
	$optimizations .= " -g3 " if ($build eq "callgrind" || $build eq "drd");
	my $strip = "strip ";
	$strip = " true " if ($build eq "debug" or $build eq "callgrind" or $build eq "drd");


	my $cppflags= " -IEngine  ";
	$cppflags .= "  -I$PsimagLite/src -I$PsimagLite $usePthreadsOrNot $floating";

	Make::make($fh,\@drivers,"DMRG++",$platform,$mpi,"$lapack $pthreadsLib",
	"$compiler $optimizations",$cppflags,$strip,"Engine/Version.h",
	"Engine/Version.h gitrev","operator");
	local *FH = $fh;
print FH<<EOF;

Engine/Version.h: gitrev
	./gitrev > Engine/Version.h

operator: dmrg
	cp dmrg operator

../doc/manual.pdf: ../doc/manual.tex
	cd ../doc; pdflatex manual.tex; pdflatex manual.tex; pdflatex manual.tex

../doc/manual.tex: ../doc/manual.ptex
	cd ../doc; find ../src -iname "*.h" -or -iname "*.cpp" | ../../PsimagLite/scripts/doc.pl manual.ptex

EOF

	close($fh);
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

