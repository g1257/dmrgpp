#!/usr/bin/perl
=pod
Copyright (c) 2009-2015, UT-Battelle, LLC
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

my ($arg) = @ARGV;

if (defined($arg) and -r "$arg" and $arg ne "Config.make") {
	my $cmd = "cp Config.make Config.make.bak";
	system($cmd);
	$cmd = "cp $arg Config.make";
	system($cmd);
}

my %dmrgMain = (name => 'dmrg', dotos => 'dmrg.o DmrgDriver1.o Provenance.o\
 RestartStruct.o FiniteLoop.o Utils.o ProgramGlobals.o');
my %dmrgDriver = (name => 'DmrgDriver1', aux => 1);
my %provenanceDriver = (name => 'Provenance', aux => 1);
my %progGlobalsDriver = (name => 'ProgramGlobals', aux => 1);
my %restartDriver = (name => 'RestartStruct', aux => 1);
my %finiteLoopDriver = (name => 'FiniteLoop', aux => 1);
my %utilsDriver = (name => 'Utils', aux => 1);
my @drivers = (\%dmrgMain,\%dmrgDriver,\%provenanceDriver,
\%progGlobalsDriver,\%restartDriver,\%finiteLoopDriver,\%utilsDriver,
"observe","toolboxdmrg");

createMakefile();

sub createMakefile
{
	unlink("Engine/Version.h");
	Make::backupMakefile();
	if (!(-r "Config.make")) {
		my $cmd = "cp Config.make.sample Config.make";
		system($cmd);
		print STDERR "$0: Executed $cmd\n";
	}

	my $fh;
	open($fh,">Makefile") or die "Cannot open Makefile for writing: $!\n";

	Make::newMake($fh,\@drivers,"DMRG++","Engine/Version.h",
	"Engine/Version.h gitrev","operator");
	local *FH = $fh;
print FH<<EOF;

gitrev: gitrev.cpp
	\$(CXX) \$(CPPFLAGS) gitrev.cpp -o gitrev

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

