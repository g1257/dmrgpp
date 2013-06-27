#!/usr/bin/perl
=pod
Copyright (c) 2009-2013, UT-Battelle, LLC
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

package Make;

sub make
{
	local *FH = shift;
	my ($drivers,$code,$platform,$mpi,$libs,$cxx,$cppflags,$strip,$additional,$additional2) = @_;
	my $allExecutables = combineAllDrivers($drivers,"");
	my $allCpps = combineAllDrivers($drivers,".cpp");

print FH<<EOF;
# DO NOT EDIT!!! Changes will be lost. Modify $0 instead
# This Makefile was written by $0
# $code by G.A.
# Platform: $platform
# MPI: $mpi

LDFLAGS =    $libs
CPPFLAGS = $cppflags
CXX = $cxx
all: $allExecutables

EOF

foreach my $what (@$drivers) {
print FH<<EOF;
$what.o: $what.cpp  Makefile $additional
	\$(CXX) \$(CPPFLAGS) -c $what.cpp

$what: $what.o
	\$(CXX) -o  $what $what.o \$(LDFLAGS)
	$strip $what
EOF
}

print FH<<EOF;

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

1;

