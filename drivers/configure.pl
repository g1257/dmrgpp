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

use lib "../../PsimagLite/scripts";
use Make;

my @drivers = ();

createMakefile();

sub createMakefile
{
	Make::backupMakefile();
	if (!(-r "Config.make")) {
		my $cmd = "cp Config.make.sample Config.make";
		system($cmd);
		print STDERR "$0: Executed $cmd\n";
	}

	my $fh;
	open($fh,">Makefile") or die "Cannot open Makefile for writing: $!\n";

	local *FH = $fh;
	my @units = qw(integrator sparseSolverTest testCRSMatrix combineContinuedFraction
	continuedFractionCollection gitrev jsonExample range kernelPolynomial
	linearPrediction options randomTest svd testLapack threads loadImbalance testIsClass
	testMemResolv1 sumDecomposition calculator closuresTest base64test);
	my $combinedUnits = combine("",\@units,".o ");
	my $combinedUnits2 = combine("./",\@units,".cpp ");

	print FH<<EOF;
include Config.make
all: @units
EOF

	foreach my $unit (@units) {
		my $doth = "../src/".ucfirst($unit).".h";
		my $tmp = (-r "$doth") ? "$doth" : "";
		print FH<<EOF;
$unit: ./$unit.cpp $tmp Makefile Makefile.dep
\t\$(CXX) \$(CPPFLAGS) -I../src -I.. -o $unit ./$unit.cpp \$(LDFLAGS)
EOF
	}

print FH<<EOF;
Makefile.dep: $combinedUnits2
\t\$(CXX) \$(CPPFLAGS) -I../src -I.. -MM  $combinedUnits2  > Makefile.dep

clean: Makefile.dep
\trm -f core* *.o *.dep *.a @units

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

