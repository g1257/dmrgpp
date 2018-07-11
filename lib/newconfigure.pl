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

use Getopt::Long qw(:config no_ignore_case);
use lib "../../PsimagLite/scripts";
use NewMake;
use PsiTag;

my @drivers = ();

my ($flavor, $config) = ("production", "");
my $usage = "USAGE: $0 [-f flavor] [-c config]\n";

GetOptions('f=s' => \$flavor,
           'c=s' => \$config) or die "$usage\n";

my @configFiles = ("../TestSuite/inputs/ConfigBase.psiTag");
push @configFiles, "../TestSuite/inputs/BasicFlavors.psiTag";
push @configFiles, $config if ($config ne "");
push @configFiles, "../TestSuite/inputs/BasicFlavors.psiTag";

createMakefile(\@configFiles, $flavor);

sub createMakefile
{
	my ($configFiles, $flavor) = @_;
	NewMake::backupMakefile();

	my $fh;
	open($fh, ">", "Makefile") or die "Cannot open Makefile for writing: $!\n";

	local *FH = $fh;
	my @units = ("MersenneTwister","Matrix","Mpi","Concurrency",
	"ProgressIndicator","MemResolv","PsimagLite","PsiBase64",
	"SpecialFunctions", "Io/TypeToH5");
	my $combinedUnits = combine("",\@units,".o ");
	my $combinedUnitsModif = $combinedUnits;
	$combinedUnitsModif =~ s/Io\///;
	my $combinedUnits2 = combine("../src/",\@units,".cpp ");
	my $configContent = NewMake::getConfigContent($configFiles, $flavor);

	print FH<<EOF;
# DO NOT EDIT!!! Changes will be lost. Use the PsiTag system to configure instead
# This Makefile was written by $0
# PsimagLite/lib

$configContent

all: libpsimaglite.a
libpsimaglite.a: Makefile $combinedUnits AinurSpirit.o AinurState.o
\tar rcs libpsimaglite.a $combinedUnitsModif AinurSpirit.o AinurState.o
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


