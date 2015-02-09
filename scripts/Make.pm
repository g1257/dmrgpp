#!/usr/bin/perl
=pod
Copyright (c) 2009-2014, UT-Battelle, LLC
All rights reserved

[PsimagLite, Version 1.0.0]

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
use File::Temp;

package Make;

sub newMake
{
	local *FH = shift;
	my ($drivers,
	    $code,
		$additional,
		$additional2,
		$additional3) = @_;
	$additional3 = "" unless defined($additional3);
	my $allExecutables = combineAllDrivers($drivers,"");
	my $allCpps = combineAllDrivers($drivers,".cpp");

print FH<<EOF;
# DO NOT EDIT!!! Changes will be lost. Modify Config.make instead
# This Makefile was written by $0
# $code by G.A.

include Config.make
CPPFLAGS += -I../../PsimagLite -I../../PsimagLite/src -IEngine
all: $allExecutables $additional3

EOF

foreach my $what (@$drivers) {
print FH<<EOF;
$what.o: $what.cpp  Makefile $additional
	\$(CXX) \$(CPPFLAGS) -c $what.cpp

$what: $what.o
	\$(CXX) -o  $what $what.o \$(LDFLAGS)
	\$(STRIP_COMMAND) $what

EOF
}

print FH<<EOF;

../../PsimagLite/lib/libpsimaglite.a:
	\$(MAKE) -f Makefile -C ../../PsimagLite/lib/

Makefile.dep: $allCpps $additional
	\$(CXX) \$(CPPFLAGS) -MM $allCpps  > Makefile.dep

clean: Makefile.dep
	rm -f core* $allExecutables *.o *.dep $additional2

include Makefile.dep
EOF
}


sub make
{
	local *FH = shift;
	my ($drivers,
	    $code,
		$platform,
		$mpi,
		$libs,
		$cxx,
		$cppflags,
		$strip,
		$additional,
		$additional2,
		$additional3) = @_;
	$additional3 = "" unless defined($additional3);
	my $allExecutables = combineAllDrivers($drivers,"");
	my $allCpps = combineAllDrivers($drivers,".cpp");

	my $gccVersion = gccVersion();
	my $normalFlags = "-Werror -Wall";
	$normalFlags .= " -Wstrict-overflow=5 " if ($gccVersion >= 4.2);
	$normalFlags .= " -frecord-gcc-switches " if ($gccVersion >= 4.3);

	my $libTarget = "";
	if ($libs=~/\-lpsimaglite/) {
		$libTarget  = " ../../PsimagLite/lib/libpsimaglite.a";
		psimagLiteLibMake($platform,$mpi,$libs,$normalFlags,$cppflags,$cxx);
	}

print FH<<EOF;
# DO NOT EDIT!!! Changes will be lost. Modify $0 instead
# This Makefile was written by $0
# $code by G.A.
# Platform: $platform
# MPI: $mpi

LDFLAGS = -L../../PsimagLite/lib   $libs
CPPFLAGS = $normalFlags $cppflags
CXX = $cxx
all: $allExecutables $additional3

EOF

foreach my $what (@$drivers) {
print FH<<EOF;
$what.o: $what.cpp  Makefile $additional $libTarget
	\$(CXX) \$(CPPFLAGS) -c $what.cpp

$what: $what.o
	\$(CXX) -o  $what $what.o \$(LDFLAGS)
	$strip $what

EOF
}

print FH<<EOF;

../../PsimagLite/lib/libpsimaglite.a:
	\$(MAKE) -f Makefile -C ../../PsimagLite/lib/

Makefile.dep: $allCpps $additional
	\$(CXX) \$(CPPFLAGS) -MM $allCpps  > Makefile.dep

clean: Makefile.dep
	rm -f core* $allExecutables *.o *.dep $additional2

include Makefile.dep
EOF
}


sub combineAllDrivers
{
	my ($drivers,$extension) = @_;
	my $buffer = "";
	foreach my $what (@$drivers) {
		my $tmp = $what.$extension." ";
		$buffer .= $tmp;
	}
	return $buffer;
}

sub findLapack
{
	my $stringToTest = "-llapack -lblas";
	my $ret = tryWith($stringToTest);
	return $stringToTest if ($ret == 0);

	$stringToTest = "/usr/lib64/liblapack.so.3 /usr/lib64/libblas.so.3";
	$ret = tryWith($stringToTest);
	return $stringToTest if ($ret == 0);

	return "/usr/lib/liblapack.so.3 /usr/lib/libblas.so.3";
}

sub tryWith
{
	my ($stringToTest) = @_;
	my $tmpfile = `mktemp XXXXXXXXXX.cpp`;
	open(FOUT,">$tmpfile") or return 1;
	print FOUT "int main() {}\n";
	close(FOUT);
	chomp($tmpfile);
	my $ret = open(PIPE,"g++ $tmpfile $stringToTest  2>&1 | ");
	if (!$ret) {
	 	system("rm -f $tmpfile");
		return 1;
	}
	my $buffer = "";
	while(<PIPE>) {
		$buffer .= $_;
	}
	close(PIPE);
	system("rm -f $tmpfile");
	$buffer =~ s/ //;
	$ret = ($buffer eq "" or $buffer eq "\n") ? 0 : 1;
	return $ret;
}

sub gccVersion
{
	my $version = 4.1;
	open(PIPE,"gcc --version |") or return $version;
	while (<PIPE>) {
		if (/^gcc +\([^\)]+\) +([^ ]+)/) {
			$version = $1;
			last;
		}
	}

	close(PIPE);

	my @temp = split(/\./,$version);
	my $n = scalar(@temp);
	my $factor = 1.0;
	$version = 0;
	for (my $i = 0; $i < $n; ++$i) {
		$version += $temp[$i]*$factor;
		$factor *= 0.1;
	}

	return $version;
}

sub psimagLiteLibMake
{
	print STDERR "$0: Make sure to compile PsimagLite/lib first\n";
}

sub backupMakefile
{
	my ($dir) = @_;
	$dir = "." unless defined($dir);
	system("cp $dir/Makefile $dir/Makefile.bak") if (-r "$dir/Makefile");
	print "Backup of $dir/Makefile in $dir/Makefile.bak\n";
}

sub findGsl
{
	my $gslDefine = " -DUSE_GSL ";
	my $gslLibs = " -lgsl -lgslcblas ";
	my $slashTmp = "/tmp";
	my @nothingFound = (" ", " ");
	return @nothingFound unless (-w $slashTmp);

	my $dir = File::Temp::tempdir(CLEANUP => 1);
	my ($fh, $filename) = File::Temp::tempfile(DIR => $dir);

	if (!$fh) {
		return @nothingFound;
	}

print $fh <<EOF;
#include "GslWrapper.h"
int main() { return 0;}
EOF
	close($fh);
	my $cppFile = $filename.".cpp";
	system("mv $filename $cppFile");
	unlink("a.out");
	system("g++ -I../../PsimagLite/src $gslDefine $cppFile  $gslLibs 2>/dev/null");
	return ($gslDefine, $gslLibs) if (-x "a.out");
	return @nothingFound;
}

1;

