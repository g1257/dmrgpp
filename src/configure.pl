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
use DmrgDriver;

my ($flavor) = @ARGV;

$flavor = procFlavor($flavor);

my %provenanceDriver = (name => 'Provenance', aux => 1);
my %progGlobalsDriver = (name => 'ProgramGlobals', aux => 1);
my %restartDriver = (name => 'RestartStruct', aux => 1);
my %finiteLoopDriver = (name => 'FiniteLoop', aux => 1);
my %utilsDriver = (name => 'Utils', aux => 1);
my %su2RelatedDriver = (name => 'Su2Related', aux => 1);
my %toolboxDriver = (name => 'toolboxdmrg',
                     dotos => 'toolboxdmrg.o ProgramGlobals.o Provenance.o Utils.o');
my $dotos = "observe.o ProgramGlobals.o Provenance.o Utils.o Su2Related.o";
my %observeDriver = (name => 'observe', dotos => $dotos);
my %kronDriver = (name => 'kronecker', dotos => 'kronecker.o ProgramGlobals.o Provenance.o');

my @drivers = (\%provenanceDriver,\%su2RelatedDriver,
\%progGlobalsDriver,\%restartDriver,\%finiteLoopDriver,\%utilsDriver,
\%observeDriver,\%toolboxDriver,\%kronDriver);

$dotos = "dmrg.o Provenance.o RestartStruct.o FiniteLoop.o Utils.o ";
$dotos .= " ProgramGlobals.o Su2Related.o";

my $templates = DmrgDriver::createTemplates();

for (my $i = 0; $i < $templates; ++$i) {
	my $name = "DmrgDriver$i";
	my %dmrgDriver = (name => $name, aux => 1);
	push @drivers,\%dmrgDriver;
	$dotos .= " $name.o ";
}

my %dmrgMain = (name => 'dmrg', dotos => $dotos);

push @drivers,\%dmrgMain;

createMakefile($flavor);

sub createMakefile
{
	my ($flavor) = @_;
	unlink("Engine/Version.h");
	Make::backupMakefile();
	Make::createConfigMake($flavor);

	my $fh;
	open($fh,">Makefile") or die "Cannot open Makefile for writing: $!\n";

	Make::newMake($fh,\@drivers,"DMRG++"," "," ","operator");
	local *FH = $fh;
print FH<<EOF;

operator: dmrg
	cp dmrg operator

../doc/manual.pdf: ../doc/manual.tex
	cd ../doc; pdflatex manual.tex; pdflatex manual.tex; pdflatex manual.tex

../doc/manual.tex: ../doc/manual.ptex
	cd ../doc; find ../src -iname "*.h" -or -iname "*.cpp" | ../../PsimagLite/scripts/doc.pl manual.ptex

EOF

	close($fh);
	print STDERR "$0: File Makefile has been written\n";
}

sub procFlavor
{
	my ($flavor) = @_;
	defined($flavor) or $flavor = "production";

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

