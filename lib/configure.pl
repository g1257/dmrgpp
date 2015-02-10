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
print FH<<EOF;
include Config.make
all: libpsimaglite.a
libpsimaglite.a: Makefile MersenneTwister.o Matrix.o Mpi.o
\tar rcs libpsimaglite.a MersenneTwister.o Matrix.o Mpi.o

MersenneTwister.o: ../src/MersenneTwister.cpp ../src/MersenneTwister.h Makefile
\t\$(CXX) \$(CPPFLAGS) -c ../src/MersenneTwister.cpp

Matrix.o: ../src/Matrix.cpp ../src/Matrix.h Makefile
\t\$(CXX) \$(CPPFLAGS) -c ../src/Matrix.cpp

Mpi.o: ../src/Mpi.cpp ../src/Mpi.h Makefile
\t\$(CXX) \$(CPPFLAGS) -c ../src/Mpi.cpp

Makefile.dep: ../src/MersenneTwister.cpp
\t\$(CXX) \$(CPPFLAGS) -MM  ../src/MersenneTwister.cpp  > Makefile.dep

clean: Makefile.dep
\trm -f core* *.o *.dep *.a

include Makefile.dep

EOF

	close($fh);
	print STDERR "File Makefile has been written\n";
}


