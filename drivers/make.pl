#!/usr/bin/perl

use strict;
use warnings;

my @drivers = ("sparseSolverTest", "testCRSMatrix", "rungeKuttaTest", "combineContinuedFraction",
"continuedFractionCollection", "gitrev", "jsonExample", "range",
"kernelPolynomial", "linearPrediction", "options", "randomTest", "svd", "testLapack", "threads");

my $lapack = findLapack();
backupMakefile();
writeMakefile();
make();

sub make
{
	system("make");
}

sub backupMakefile
{
	system("cp Makefile Makefile.bak") if (-r "Makefile");
	print "Backup of Makefile in Makefile.bak\n";
}

sub writeMakefile
{
	my $allExecutables = combineAllDrivers("");
	my $allCpps = combineAllDrivers(".cpp");

open(FILE,">Makefile") or die "Cannot open Makefile for writing: $!\n";
print FILE<<EOF;
# DO NOT EDIT!!! Changes will be lost. Modify configure.pl instead
# This Makefile was written by configure.pl
# PsimagLite by G.A.
# Platform: Linux
# MPI: 0

LDFLAGS =      $lapack    -lm  -lpthread 
CPPFLAGS = -Werror -Wall -I../  -I../src
CXX = g++ -O3 -DNDEBUG

all: $allExecutables
EOF

foreach my $what (@drivers) {
print FILE<<EOF;
$what.o: $what.cpp  Makefile
	\$(CXX) \$(CPPFLAGS) -c $what.cpp  

$what: $what.o
	\$(CXX) -o  $what $what.o \$(LDFLAGS)

EOF
}

print FILE<<EOF;
Makefile.dep: $allCpps
	\$(CXX) \$(CPPFLAGS) -MM $allCpps  > Makefile.dep

clean:
	rm -f core* $allExecutables *.o *.ii *.tt

include Makefile.dep
######## End of Makefile ########

EOF

close(FILE);
print "Done writing Makefile\n";
}

sub combineAllDrivers
{
	my ($extension) = @_;
	my $buffer = "";
	foreach my $what (@drivers) {
		my $tmp = $what.$extension." ";
		$buffer .= $tmp;
	}
	return $buffer;
}

sub findLapack
{
	my $stringToTest = "-llapack -lblas";
	my $ret = tryWith($stringToTest);
	return $stringToTest if ($ret == 0);

	$stringToTest = "/usr/lib64/liblapack.so.3 /usr/lib64/libblas.so.3";
	$ret = tryWith($stringToTest);
	return $stringToTest if ($ret == 0);

	return "/usr/lib/liblapack.so.3 /usr/lib/libblas.so.3";
}

sub tryWith
{
	my ($stringToTest) = @_;
	my $tmpfile = `mktemp`;
	open(FOUT,">$tmpfile") or return 1;
	print FOUT "int main() {}\n";
	close(FOUT);
	my $ret = open(PIPE,"g++ test.cpp $stringToTest  2>&1 | ");
	if (!$ret) {
	 	system("rm -f $tmpfile");	
		return 1;
	}
	my $buffer = "";
	while(<PIPE>) {
		$buffer .= $_;
	}
	close(PIPE);
	system("rm -f $tmpfile");
	$buffer =~ s/ //;
	$ret = ($buffer eq "" or $buffer eq "\n") ? 0 : 1;
	return $ret;
}


