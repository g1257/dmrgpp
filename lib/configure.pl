#!/usr/bin/perl
=pod
Copyright (c) 2009-2015, UT-Battelle, LLC
All rights reserved

[PsimagLite, Version 1.0]

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

my @drivers = ();

my ($flavor) = @ARGV;

$flavor = procFlavor($flavor);
createMakefile($flavor);

sub createMakefile
{
	my ($flavor) = @_;
	Make::backupMakefile();
	Make::createConfigMake($flavor);

	if (!(-r "Config.make")) {
		my $cmd = "cp Config.make.sample Config.make";
		system($cmd);
		print STDERR "$0: Executed $cmd\n";
	}

	my $fh;
	open($fh,">Makefile") or die "Cannot open Makefile for writing: $!\n";

	local *FH = $fh;
	my @units = ("MersenneTwister","Matrix","Mpi","Concurrency",
	"ProgressIndicator","MemResolv","PsimagLite","PsiBase64",
	"SpecialFunctions");
	my $combinedUnits = combine("",\@units,".o ");
	my $combinedUnits2 = combine("../src/",\@units,".cpp ");

	print FH<<EOF;
include Config.make
all: libpsimaglite.a
libpsimaglite.a: Makefile $combinedUnits AinurSpirit.o AinurState.o
\tar rcs libpsimaglite.a $combinedUnits AinurSpirit.o AinurState.o
EOF

	foreach my $unit (@units) {
		print FH<<EOF;
$unit.o: ../src/$unit.cpp ../src/$unit.h Makefile
\t\$(CXX) \$(CPPFLAGS) -c ../src/$unit.cpp
EOF
	}

print FH<<EOF;
AinurSpirit.o: ../src/Ainur/AinurSpirit.cpp ../src/Ainur/AinurSpirit.h Makefile ../src/Ainur/AinurSpirit1.cpp
\t\$(CXX) \$(CPPFLAGS) -c ../src/Ainur/AinurSpirit.cpp

AinurState.o: ../src/Ainur/AinurState.cpp ../src/Ainur/AinurState.h Makefile 
\t\$(CXX) \$(CPPFLAGS) -c ../src/Ainur/AinurState.cpp

Makefile.dep: $combinedUnits2 ../src/Ainur/AinurSpirit.cpp ../src/Ainur/AinurState.cpp
\t\$(CXX) \$(CPPFLAGS) -MM  $combinedUnits2  > Makefile.dep

clean: Makefile.dep
\trm -f core* *.o *.dep *.a

include Makefile.dep

EOF

	close($fh);
	print STDERR "File Makefile has been written\n";
}

sub combine
{
	my ($pre,$a,$post) = @_;
	my $n = scalar(@$a);
	my $buffer = "";
	for (my $i = 0; $i < $n; ++$i) {
		$buffer .= $pre.$a->[$i].$post;
	}

	return $buffer;
}

sub procFlavor
{
	my ($flavor) = @_;
	if (!defined($flavor)) {
		$flavor = "production";
		print STDERR "$0: No flavor given, assuming production\n";
		print STDERR "\t say $0 help for a list of options\n";
	}

	my $hasPath = ($flavor =~ /^\.\./ or $flavor =~ /^\//);
	return $flavor if ($hasPath);

	if ($flavor eq "help") {
		print "USAGE: $0 [production | debug | callgrind";
		print " | helgrind | drd]\n";
		exit(0);
	}

	my $dir = "../TestSuite/inputs";
	if ($flavor eq "production") {
		$flavor = "Config.make";
	} elsif ($flavor eq "debug") {
		$flavor = "ConfigDebug.make";
	} elsif ($flavor eq "callgrind") {
		$flavor = "ConfigCallgrind.make";
	} elsif ($flavor eq "helgrind" or $flavor eq "drd") {
		$flavor = "ConfigHelgrind.make";
	} else {
		return $flavor;
	}

	return "$dir/$flavor";
}

