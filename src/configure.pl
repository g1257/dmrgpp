#!/usr/bin/perl
=pod
Copyright (c) 2020, UT-Battelle, LLC
All rights reserved

[Dmft, Version 0.]

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
use PsiTag;

my ($flavor, $lto) = (NewMake::noFlavor(), 0);
my $usage = "USAGE: $0 [-f flavor] [-lto] [-c config]\n";
my $config;

GetOptions('f=s' => \$flavor,
           'lto' => \$lto,
           'c=s' => \$config) or die "$usage\n";

my $gccdash = "";
if ($lto == 1) {
	$gccdash = "gcc-";
	$lto = "-flto";
} else {
	$lto = "";
}

my $basicConfig = "ConfigBase.psiTag";
my @configFiles = NewMake::configFilesList($basicConfig, $config);


my $libs = "../../dmrgpp/src/libdmrgpp.a ../../dmrgpp/src/KronUtil/libkronutil.a";
my $dotos = "cincuenta.o $libs";

my %cincuentaMain = (name => 'cincuenta', dotos => "$dotos");
my %provenance = (name => 'Provenance', aux => 1);

my @drivers = (\%cincuentaMain, \%provenance);

my %args;
$args{"CPPFLAGS"} = $lto;
$args{"LDFLAGS"} = $lto;
$args{"flavor"} = $flavor;
$args{"code"} = "Dmft";
$args{"configFiles"} = \@configFiles;
$args{"additional3"} = "GitRevision.h";

system("./createGitRevision.pl GitRevision.h");

createMakefile(\@drivers, \%args);

sub createMakefile
{
	my ($drivers, $args) = @_;
	unlink("Makefile.dep");
	NewMake::backupMakefile();

	my $fh;
	open($fh, ">", "Makefile") or die "Cannot open Makefile for writing: $!\n";
	local *FH = $fh;
print FH<<EOF;

.PHONY: GitRevision.h

GitRevision.h:
	./createGitRevision.pl GitRevision.h
EOF

	NewMake::main($fh, $args, $drivers);
	close($fh);

	print STDERR "$0: File Makefile has been written\n";
}

