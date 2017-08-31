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

my ($flavor, $generateSources) = @ARGV;
$flavor = procFlavor($flavor);
my $gccdash = "";
my $lto = "";
if (defined($generateSources)) {
	if ($generateSources eq "lto") {
		$gccdash = "gcc-";
		$lto = "-flto";
		$generateSources = 0;
	}
} else {
	$generateSources = 0;
}

system("cd KronUtil; perl configure.pl $gccdash");

my %provenanceDriver = (name => 'Provenance', aux => 1);
my %progGlobalsDriver = (name => 'ProgramGlobals', aux => 1);
my %restartDriver = (name => 'RestartStruct', aux => 1);
my %finiteLoopDriver = (name => 'FiniteLoop', aux => 1);
my %utilsDriver = (name => 'Utils', aux => 1);
my %su2RelatedDriver = (name => 'Su2Related', aux => 1);
my %toolboxDriver = (name => 'toolboxdmrg',
                     dotos => 'toolboxdmrg.o ProgramGlobals.o Provenance.o Utils.o');
my $dotos = "observe.o ProgramGlobals.o Provenance.o Utils.o Su2Related.o";
$dotos .= " ObserveDriver0.o ObserveDriver1.o ObserveDriver2.o ";
my %observeDriver = (name => 'observe', dotos => $dotos);
my %kronDriver = (name => 'kronecker', dotos => 'kronecker.o ProgramGlobals.o Provenance.o');

my %observeDriver0 = (name => 'ObserveDriver0', aux => 1);
my %observeDriver1 = (name => 'ObserveDriver1', aux => 1);
my %observeDriver2 = (name => 'ObserveDriver2', aux => 1);

my @drivers = (\%provenanceDriver,\%su2RelatedDriver,
\%progGlobalsDriver,\%restartDriver,\%finiteLoopDriver,\%utilsDriver,
\%observeDriver,\%toolboxDriver,\%kronDriver,
\%observeDriver0,\%observeDriver1,\%observeDriver2);

$dotos = "dmrg.o Provenance.o RestartStruct.o FiniteLoop.o Utils.o ";
$dotos .= " ProgramGlobals.o Su2Related.o";

my $templates = DmrgDriver::createTemplates($generateSources);

for (my $i = 0; $i < $templates; ++$i) {
	my $name = "DmrgDriver$i";
	my %dmrgDriver = (name => $name, aux => 1);
	push @drivers,\%dmrgDriver;
	$dotos .= " $name.o ";
}

my %dmrgMain = (name => 'dmrg', dotos => "$dotos", libs => "kronutil");

push @drivers,\%dmrgMain;

createMakefile($flavor, $lto);

sub createMakefile
{
	my ($flavor, $lto) = @_;
	Make::backupMakefile();
	my %args;
	$args{"CPPFLAGS"} = $lto;
	$args{"LDFLAGS"} = $lto;
	Make::createConfigMake($flavor, \%args);

	my $fh;
	open($fh, ">", "Makefile") or die "Cannot open Makefile for writing: $!\n";

	my %additionals;
	$additionals{"code"} = "DMRG++";
	$additionals{"additional3"} = "operator";
	Make::newMake($fh,\@drivers,\%additionals);
	local *FH = $fh;
print FH<<EOF;

operator: dmrg
	cp dmrg operator

libkronutil.a:
	\$(MAKE) -C KronUtil

../doc/manual.pdf: ../doc/manual.tex
	cd ../doc; pdflatex manual.tex; pdflatex manual.tex; pdflatex manual.tex

../doc/manual.tex: ../doc/manual.ptex
	cd ../doc; find ../src -iname "*.h" -or -iname "*.cpp" | ../../PsimagLite/scripts/doc.pl manual.ptex

clean::
	\$(MAKE) -C KronUtil clean
EOF

	close($fh);
	print STDERR "$0: File Makefile has been written\n";
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

