#!/usr/bin/perl
=pod
Copyright (c) 2009-2017, UT-Battelle, LLC
All rights reserved

[PsimagLite, Version 1.]

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

my $flavor = NewMake::noFlavor();
my $usage = "USAGE: $0 [-f flavor]\n";
my $config;

GetOptions('f=s' => \$flavor,
           'c=s' => \$config) or die "$usage\n";

my @configFiles = ("../../dmrgpp/TestSuite/inputs/ConfigBase.psiTag");
push @configFiles, $config if (defined($config));

createMakefile(\@configFiles, $flavor);

sub createMakefile
{
	my ($configFiles, $flavor) = @_;

	my @drivers = qw(integrator sparseSolverTest testCRSMatrix combineContinuedFraction
	continuedFractionCollection range kernelPolynomial fit
	linearPrediction options randomTest svd testLapack threads loadImbalance testIsClass
	testMemResolv1 sumDecomposition calculator closuresTest base64test checkRunId
	testLanczos testExcitedLanczos testLanczosMatrixInFile nested testIoNg testIoNgBoolean
	affinityTest testPredicate isBlasThreaded);

	my %args;
	$args{"code"} = "PsimagLite/drivers";
	$args{"configFiles"} = $configFiles;
	$args{"flavor"} = $flavor;

	NewMake::backupMakefile();
	my $fh;
	open($fh, ">", "Makefile") or die "Cannot open Makefile for writing: $!\n";

	NewMake::main($fh, \%args, \@drivers);
}


