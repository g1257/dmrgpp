#!/usr/bin/perl

use strict;
use warnings;

use lib '../scripts';
use Make;

my @drivers = ("integrator","sparseSolverTest", "testCRSMatrix", "rungeKuttaTest", "combineContinuedFraction",
"continuedFractionCollection", "gitrev", "jsonExample", "range",
"kernelPolynomial", "linearPrediction", "options", "randomTest", "svd", "testLapack", "threads",
"testIsClass","testMemResolv1","sumDecomposition","calculator","closuresTest");

my $lapack = Make::findLapack();
Make::backupMakefile();
writeMakefile();
make();

sub make
{
	system("make");
}

sub writeMakefile
{
	open(my $fh,">Makefile") or die "Cannot open Makefile for writing: $!\n";

	my $libs = " -lm  -lpthread -lpsimaglite $lapack";
	my $cxx = "g++ -O3 -DNDEBUG";
	my $cppflags = " -I../  -I../src";
	Make::make($fh,\@drivers,"PsimagLite","Linux",0,$libs,$cxx,$cppflags,"true"," "," ");

	close($fh);
	print "$0: Done writing Makefile\n";
}
