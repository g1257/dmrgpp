#!/usr/bin/perl
=pod
Copyright (c) 2009-2018, UT-Battelle, LLC
All rights reserved

[DMRG++, Version 5.]

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
use lib ".";
use DmrgDriver;
use PsiTag;

my ($flavor, $generateSources, $su2enabled, $lto) = (NewMake::noFlavor() , 0, 0, 0);
my $usage = "USAGE: $0 [-f flavor] [-s] [-su2] [-lto] [-c config]\n";
my $config;

GetOptions('f=s' => \$flavor,
           's' => \$generateSources,
           'su2' => \$su2enabled,
           'lto' => \$lto,
           'c=s' => \$config) or die "$usage\n";

my $gccdash = "";
if ($lto == 1) {
	$gccdash = "gcc-";
	$lto = "-flto";
} else {
	$lto = "";
}

my $basicConfig = "../../dmrgpp/TestSuite/inputs/ConfigBase.psiTag";
my @configFiles = NewMake::configFilesList($basicConfig, $config);

system("cd KronUtil; perl configure.pl \"@configFiles\" $flavor $gccdash");

my %provenanceDriver = (name => 'Provenance', aux => 1);
my %progGlobalsDriver = (name => 'ProgramGlobals', aux => 1);
my %restartDriver = (name => 'RestartStruct', aux => 1);
my %finiteLoopDriver = (name => 'FiniteLoop', aux => 1);
my %utilsDriver = (name => 'Utils', aux => 1);
my %qnDriver = (name => 'Qn', aux => 1);
my %su2RelatedDriver = (name => 'Su2Related', aux => 1);
my %toolboxDriver = (name => 'toolboxdmrg',
                     dotos => 'toolboxdmrg.o ProgramGlobals.o Provenance.o Utils.o Qn.o');
my $dotos = "observe.o ProgramGlobals.o Provenance.o Utils.o Su2Related.o Qn.o ";
$dotos .= " ObserveDriver0.o ObserveDriver1.o ObserveDriver2.o ";
my %observeDriver = (name => 'observe', dotos => $dotos);

my %observeDriver0 = (name => 'ObserveDriver0', aux => 1);
my %observeDriver1 = (name => 'ObserveDriver1', aux => 1);
my %observeDriver2 = (name => 'ObserveDriver2', aux => 1);

my @drivers = (\%provenanceDriver,\%su2RelatedDriver,
\%progGlobalsDriver,\%restartDriver,\%finiteLoopDriver,\%utilsDriver,
\%qnDriver, \%observeDriver,\%toolboxDriver,
\%observeDriver0,\%observeDriver1,\%observeDriver2);

$dotos = "dmrg.o Provenance.o RestartStruct.o FiniteLoop.o Utils.o Qn.o ";
$dotos .= " ProgramGlobals.o Su2Related.o";

my @su2files = DmrgDriver::createTemplates($generateSources);
my $templates = scalar(@su2files);

for (my $i = 0; $i < $templates; ++$i) {
	next if (!$su2enabled && $su2files[$i]);
	my $name = "DmrgDriver$i";
	my %dmrgDriver = (name => $name, aux => 1);
	push @drivers,\%dmrgDriver;
	$dotos .= " $name.o ";
}

my %dmrgMain = (name => 'dmrg', dotos => "$dotos", libs => "kronutil");

push @drivers,\%dmrgMain;

my $su2flags = ($su2enabled) ? " -DENABLE_SU2 " : "";

my %args;
$args{"CPPFLAGS"} = $lto.$su2flags;
$args{"LDFLAGS"} = $lto;
$args{"flavor"} = $flavor;
$args{"code"} = "DMRG++";
$args{"configFiles"} = \@configFiles;
$args{"additional3"} = "GitRevision.h";
$args{"additional4"} = $args{"additional3"};

system("./createGitRevision.pl GitRevision.h");

createMakefile(\@drivers, \%args);

sub createMakefile
{
	my ($drivers, $args) = @_;
	unlink("Makefile.dep");
	NewMake::backupMakefile();
	$args->{"additional3"} = "operator";

	my $fh;
	open($fh, ">", "Makefile") or die "Cannot open Makefile for writing: $!\n";

	NewMake::main($fh, $args, $drivers);
	local *FH = $fh;
print FH<<EOF;

.PHONY: GitRevision.h

GitRevision.h:
	./createGitRevision.pl GitRevision.h

operator: dmrg
	cp dmrg operator

testQn: testQn.o Qn.o
	\$(CXX) Qn.o testQn.o \$(LDFLAGS) -o testQn

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

