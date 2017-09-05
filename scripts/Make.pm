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
	my ($drivers, $additionals) = @_;
	my %a = %$additionals;
	my $additional = $a{"additional"};
	my $additional2 = $a{"additional2"};
	my $additional3 = $a{"additional3"};
	my $path = $a{"path"};
	my $code = $a{"code"};
	$additional = " " unless defined($additional);
	$additional2 = " " unless defined($additional2);
	$additional3 = "" unless defined($additional3);
	$path = " " unless defined($path);

	my $allExecutables = combineAllDrivers($drivers,"");
	my $allCpps = combineAllDrivers($drivers,".cpp");

print FH<<EOF;
# DO NOT EDIT!!! Changes will be lost. Modify Config.make instead
# This Makefile was written by $0
# $code

include ${path}Config.make
CPPFLAGS += -I$path../../PsimagLite -I$path../../PsimagLite/src -I${path}Engine
all: $allExecutables $additional3

EOF

foreach my $ptr (@$drivers) {
	my $refptr = ref($ptr);
	my $oldmode = ($refptr eq "");
	my $what = ($oldmode) ? $ptr : $ptr->{"name"};
	my $aux = ($oldmode) ? 0 : $ptr->{"aux"};
	$aux = 0 if (!defined($aux));
	my $dotos = ($oldmode) ? "$what.o" : $ptr->{"dotos"};
	$dotos = "$what.o" if (!defined($dotos));

	print FH<<EOF;
$what.o: $what.cpp  Makefile $additional ${path}Config.make
	\$(CXX) \$(CPPFLAGS) -c $what.cpp

EOF

	if (!$aux) {
		# FIXME: Support many libs separated by commas here
		my $libs = ($oldmode) ? "" : $ptr->{"libs"};
		my $libs1 = "";
		my $libs2 = "";
		if (defined($libs) and $libs ne "") {
			$libs1 = "lib$libs.a";
			$libs2 = "-l$libs";
		}

		print FH<<EOF;
$what: $dotos $libs1
	\$(CXX) -o  $what $dotos \$(LDFLAGS) $libs2 \$(CPPFLAGS)
	\$(STRIP_COMMAND) $what

EOF
	}
}

print FH<<EOF;

$path../../PsimagLite/lib/libpsimaglite.a:
	\$(MAKE) -f Makefile -C $path../../PsimagLite/lib/

Makefile.dep: $allCpps $additional
	\$(CXX) \$(CPPFLAGS) -MM $allCpps  > Makefile.dep

clean::
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

	my $libTarget = "";
	if ($libs=~/\-lpsimaglite/) {
		$libTarget  = " ../../PsimagLite/lib/libpsimaglite.a";
		psimagLiteLibMake();
	}

print FH<<EOF;
# DO NOT EDIT!!! Changes will be lost. Modify $0 instead
# This Makefile was written by $0
# $code by G.A.
# Platform: $platform
# MPI: $mpi

LDFLAGS = -L../../PsimagLite/lib   $libs
CPPFLAGS = $cppflags
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
	foreach my $ptr (@$drivers) {
		my $refptr = ref($ptr);
		my $oldmode = ($refptr eq "");
		my $what = ($oldmode) ? $ptr : $ptr->{"name"};
		my $aux = ($oldmode) ? 0 : $ptr->{"aux"};
		defined($aux) or $aux = 0;
		next if ($aux and $extension eq "");
		my $tmp = $what.$extension." ";
		$buffer .= $tmp;
	}

	return $buffer;
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
	print STDERR "$0: Backup of $dir/Makefile in $dir/Makefile.bak\n";
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

sub createConfigMake
{
	my ($flavor, $additionals) = @_;
	defined($flavor) or $flavor = "../TestSuite/inputs/Config.make";
	my $cmd = "cp ../TestSuite/inputs/ConfigBase.make Config.make.new";
	system($cmd);
	print STDERR "$0: Executed $cmd\n";
	if (!(-r "$flavor")) {
		print STDERR "$0: WARNING: I didn't find $flavor, ignoring...\n";
		return;
	}

	open(FILE, "<", $flavor) or return;
	if (!open(FOUT, ">>", "Config.make.new")) {
		close(FILE);
		return;
	}

	my $blasLapack = " -llapack -lblas ";
	while (<FILE>) {
		next if (/ConfigBase\.make/);
		print FOUT;
	}

	close(FILE);

	if (defined($additionals)) {
		my $cppflags = $additionals->{"CPPFLAGS"};
		my $ldflags = $additionals->{"LDFLAGS"};
		print FOUT "CPPFLAGS += $cppflags\n" if ($cppflags ne "");
		print FOUT "LDFLAGS += $ldflags\n" if ($ldflags ne "");
	}

	close(FOUT);

	my $overwrite = 1;
	if (-r "Config.make") {
		$overwrite = 0;
		print "$0: Config.make exists. Do you want to overwrite it?\n";
		print "$0: Your answer: (n/Y): ";
		$_ = <STDIN>;
		chomp;
		$overwrite = ($_ eq "Y");
		if (!$overwrite) {
			print STDERR "$0: Please consider comparing your Config.make\n";
			print STDERR "\t to the one I've written to Config.make.new\n";
			return;
		}
	}

	$cmd = "mv Config.make.new Config.make";
	system($cmd);
	print STDERR "$0: Written Config.make\n";
}

1;

