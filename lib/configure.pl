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
	my @units = ("MersenneTwister","Matrix","Mpi","ApplicationInfo","Concurrency",
	"ProgressIndicator","Tokenizer");
	my $combinedUnits = combine("",\@units,".o ");
	my $combinedUnits2 = combine("../src/",\@units,".cpp ");

	print FH<<EOF;
include Config.make
all: libpsimaglite.a
libpsimaglite.a: Makefile $combinedUnits
\tar rcs libpsimaglite.a $combinedUnits
EOF

	foreach my $unit (@units) {
		print FH<<EOF;
$unit.o: ../src/$unit.cpp ../src/$unit.h Makefile
\t\$(CXX) \$(CPPFLAGS) -c ../src/$unit.cpp
EOF
	}

print FH<<EOF;
Makefile.dep: $combinedUnits2
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

