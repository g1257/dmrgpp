#!/usr/bin/perl -w

=pod
// BEGIN LICENSE BLOCK
/*
Copyright © 2008 , UT-Battelle, LLC
All rights reserved

[DMRG++, Version 2.0.0]
[by G.A., Oak Ridge National Laboratory]

UT Battelle Open Source Software License 11242008
see file LICENSE for more details
*********************************************************
DISCLAIMER

THE SOFTWARE IS SUPPLIED BY THE COPYRIGHT HOLDERS AND
CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED
WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A
PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE
COPYRIGHT OWNER, CONTRIBUTORS, UNITED STATES GOVERNMENT,
OR THE UNITED STATES DEPARTMENT OF ENERGY BE LIABLE FOR
ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE
OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH
DAMAGE.

NEITHER THE UNITED STATES GOVERNMENT, NOR THE UNITED
STATES DEPARTMENT OF ENERGY, NOR THE COPYRIGHT OWNER, NOR
ANY OF THEIR EMPLOYEES, REPRESENTS THAT THE USE OF ANY
INFORMATION, DATA, APPARATUS, PRODUCT, OR PROCESS
DISCLOSED WOULD NOT INFRINGE PRIVATELY OWNED RIGHTS.

*********************************************************

// END LICENSE BLOCK

"Program testing can be a very effective way to show the presence of bugs, 
but it is hopelessly inadequate for showing their absence." -- Edsger Dijkstra

#TAGSTART DO NOT REMOVE THIS TAG
List of Standard Tests for DMRG++ version 2.0.0
(See Notation in the doc directory)
Tests *with* SU(2) symmetry:
0) Hubbard Model One Orbital (HuStd-1orb) on a chain (CubicStd1d) for U=1 with 8+8 sites. 
	INF(100)+7(100)-7(100)-7(100)+7(100)
1) Hubbard Model One Orbital (HuStd-1orb) on a chain (CubicStd1d) for U=2 with 8+8 sites
	INF(60)+7(100)-7(100)-7(100)+7(100)
2) Hubbard Model One Orbital (HuStd-1orb) on a chain (CubicStd1d) for U=2 with 8+8 sites
        INF(100)+7(200)-7(200)-7(200)+7(200) [to calculate correlations]
3) Hubbard Model One Orbital (HuStd-1orb) on a ladder (CubicStd2d) for U=0 with 8+8 sites
	INF(100)+7(200)-7(200)-7(200)+7(200) [to calculate correlations]
4) Hubbard Model One Orbital (HuStd-1orb) on a chain (CubicStd1d) for U=1 with 8+8 sites.
        INF(60)+7(100)-7(100)-7(100)+7(100) test for TimeStepTargetting "creation op"
5) Hubbard Model One Orbital (HuStd-1orb) on a chain (CubicStd1d) for U=1 with 8+8 sites.
        INF(60)+7(100)-7(100)-7(100)+7(100) test for TimeStepTargetting "identity op"
6) Hubbard Model One Orbital (HuStd-1orb) on a chain (CubicStd1d) for U=1 V=-0.5 with 30+30 sites.
        INF(226)+29(226)-29(226)-29(226)+29(226) with angularMomentum j=5 (as in the literature)
7) Like test 4 but with many times
8) Like test 5 but with many times
20) Heisenberg Model Spin 1/2 (HeStd-F12) on a chain (CubicStd1d) for J=1 with 16+16 sites
	INF(60)+7(100)-7(100)-7(100)+7(100)
21) Heisenberg Model Spin 1/2 (HeStd-F12) on a chain (CubicStd1d) for J=2.5 with 8+8 sites
	INF(100)+7(200)-7(200)-7(200)+7(200)
22) Heisenberg Model Spin 1/2 (HeStd-F12) on a chain (CubicStd1d) for J=1.0 with 8+8 sites
	INF(100)+7(200)-7(200)-7(200)+7(200)+7(200)+1(200) [to calculate correlations]
23) Heisenberg Model Spin 1/2 (HeStd-F12) on a chain (CubicStd1d) for J=1.0 with 8+8 sites
        checkpointA
24) Heisenberg Model Spin 1/2 (HeStd-F12) on a chain (CubicStd1d) for J=1.0 with 8+8 sites
        checkpointB
25)  Heisenberg Model Spin 1/2 (HeStd-F12) on a chain (CubicStd1d) for J=2.5 with 8+8 sites
        INF(100)+7(200)-7(200)-7(200)+7(200) To check the WFT
40) Fe-based Superconductors model (HuFeAS-2orb) on a ladder (LadderFeAs) with U=0 J=0 with 4+4 sites
	 INF(60)+7(100)-7(100)-7(100)+7(100)
41) Fe-based Superconductors model (HuFeAS-2orb) on a ladder (LadderFeAs) with U=1 J=1 with 4+4 sites
	INF(60)+7(100)-7(100)-7(100)+7(100)
60) t-J one-Orbital (tjOneOrbital) on a chain with t=1 J=1 with 16+16 sites
	 INF(60)+7(100)-7(100)-7(100)+7(100)

Note 1: Tests *without* SU(2) symmetry: Add 100 to the above test.
	For example test 102) is test 2) with SU(2) symmetry *disabled*.

Note 2: This test suite does not support pthreaded tests yet.
#TAGEND DO NOT REMOVE THIS TAG
=cut

use strict;
use Getopt::Long;
my $compareTool = "kompare"; # Not needed, just in case you have it
my $executable;
my $testNumber;
my $all;
my @GlobalHashTable;
my ($generateFlag,$buildFlag,$runFlag,$cmpFlag) = (1,1,1,1);
GetOptions ("exe=s" => \$executable,"n=i" => \$testNumber,"all"  => \$all,
	"g=i" => \$generateFlag,"b=i" => \$buildFlag, "r=i" => \$runFlag,"c=i" => \$cmpFlag);

$all=0 if (!defined($all));

$testNumber = askWhatTest() if (!defined($testNumber) and $all==0);

if ($all) { # Run all tests 
	runAllTests(0);
} elsif ($testNumber<0) { # Run all tests starting from -$testNumber
	my $start = -$testNumber;
	runAllTests($start);
} else { # Run test number $testNumber
	testSuite($testNumber,$executable);
}

sub runAllTests
{
	my ($start)=@_;
	my $available = getAvailableTests(0);
	my @temp = split(/ /,$available);
	my $execSaved;
	for (my $i=0;$i<$#temp+1;$i++) {
		next if ($temp[$i] eq "");
		next if ($temp[$i]<$start);
		my $exec;
		$exec = $execSaved if ($i>0 && $temp[$i]!=100 && thisSpecEqualToPrevSpec($temp[$i],$temp[$i-1]));
		
		$execSaved = testSuite($temp[$i],$exec);
	}
	print "All tests done successfully.\n";
}

sub testSuite
{
	my ($tn,$executable) = @_;
	
	return "nothing" if (skipMe($tn));
	
	$executable = createExecutable($tn) if (!defined($executable));

	runTest($tn,$executable) if ($runFlag);

	if ($cmpFlag) {
		validateEnergy($tn);
		validateProfile($tn);
	}
	return $executable;
}

sub skipMe
{
	my ($tn)=@_;
	open(PIPE,"inputs/input$tn.inp") or die "Cannot open inputs/input$tn.inp: $!\n";
	$_=<PIPE>;
	close(PIPE);
	return 1 if (/^skip/i);
	return 0;
}

sub runTest
{
	my ($tn,$executable) = @_;
	print "Please Wait while the test is run...\n";
	unlink("data$tn.txt"); 
	my $rCode = system("$executable inputs/input$tn.inp >& /dev/null ");
	die "$0: Test cannot be run to completion\n" if ($rCode != 0);
	print "The run has been completed.\n";
}

sub validateEnergy
{
	my ($testNumber) = @_;
	system("grep Energy data$testNumber.txt > resultsOfTest$testNumber.txt");
	system("diff oracles/validatedEnergy$testNumber.txt resultsOfTest$testNumber.txt > diffOfTest$testNumber.txt ");
	my $rCode = system("$compareTool validated$testNumber.txt resultsOfTest$testNumber.txt >& /dev/null");
	if ($rCode != 0) {	
		print "See file diffOfTest$testNumber.txt\n";
	}
}

sub validateProfile
{
	my ($testNumber) = @_;
	system("gprof ../src/dmrg > profileOfTest$testNumber.txt");
	system("diff profileOfTest$testNumber.txt oracles/validateProfile$testNumber.txt > diffOfProfile$testNumber.txt");
}


sub validateChargeCorrelations
{
	my ($testNumber) = @_;
	system("../src/observe input$testNumber.inp data$testNumber.txt  > observeData$testNumber.txt");
	system("diff oracles/validateObserve$testNumber.txt observeData$testNumber.txt");
	
}

sub askWhatTest
{
	my $available = getAvailableTests(1);
	print "Please type the number of the test you want to run\n";
	print "Available: $available\n";
	print "Default: 0 (press ENTER) ";
	$_ = <STDIN>;
	chomp;
	s/ //g;
	$_ = "0" if ($_ eq "");
	print "OK, I'll run test number $_\n";
	return $_;
}

sub createExecutableOld
{
	my ($testNumber) = @_;
	my $tn = $testNumber;
	$tn -= 100 if ($tn>=100);
	print "Trying to create executable for test $testNumber...\n";
	chdir("../src");
	my $specFile = getSpecFile($tn);
	my $rCode = 0;
	if ($generateFlag) {
		$rCode = system("perl configure.pl < $specFile >& /dev/null"); 
		die "Problem creating executable (configure.pl failed)\n" if ($rCode != 0);
	}
	if ($buildFlag) {
		$rCode = system("make");
		die "Problem creating executable (make failed)\n" if ($rCode != 0);
	}
	chdir ("../TestSuite");
	return "nice -10 ../src/dmrg";
}

sub createExecutable
{
	my ($testNumber) = @_;
	my ($exeIndex,$h);
	($exeIndex,$h) = getExecHashAndIndex($testNumber);
	
	if ($exeIndex<0) {
		print "Trying to create executable for test $testNumber...\n";
		chdir("../src");
		my $specFile = getSpecFile($testNumber);
		my $rCode = 0;
	
		if ($generateFlag) {
			$rCode = system("perl configure.pl < $specFile >& /dev/null"); 
			die "Problem creating executable (configure.pl failed)\n" if ($rCode != 0);
		}
		if ($buildFlag) {
			$rCode = system("make");
			die "Problem creating executable (make failed)\n" if ($rCode != 0);
			system("cp dmrg dmrg-$h");
		}
		chdir ("../TestSuite");
	}
	
	return "nice -10 ../src/dmrg-$h";
}


sub getExecHashAndIndex
{
		my ($tn)=@_;
		my $specFile = getSpecFile($tn);
		open(PIPE,"md5sum $specFile |") or die "Cannot open pipe: $!\n";
		$_ = <PIPE>;
		my @temp = split;
		close(PIPE);
		my $i = getExeIndex($temp[0]);
		return ($i,$temp[0]);
}

sub getExeIndex
{
	my ($h) = @_;
	my $n = $#GlobalHashTable+1;
	for (my $i=0;$i<$n;$i++) {
		return $i if ($GlobalHashTable[$i] eq $h);
	}
	$GlobalHashTable[$n] = $h;
	return -1;
}


sub getSpecFile
{
	my ($tn) = @_;
	my $t = $tn;
	$t = $tn -100 if ($tn>=100);
	return "../TestSuite/inputs/model$t.spec";
}

sub thisSpecEqualToPrevSpec
{
	my ($tnNew,$tnOld)=@_;
	my $fNew = getSpecFile($tnNew);
	my $fOld = getSpecFile($tnOld);
	return filesAreEqual($fNew,$fOld);
}

sub filesAreEqual
{
	my ($f1,$f2) = @_;
	#print "Comparing $f1 with $f2\n";
	open(PIPE,"diff $f1 $f2 | wc -l | ") or return 0;
	$_ = <PIPE>;
	close(PIPE);
	return 1 if (/^0$/);
	#die "Files not equal $f1 $f2 $_\n";
	return 0;
}

sub getAvailableTests
{
	my ($needsPrinting)=@_; # print all available tests?
	
	open(MYSELF,$0) or die "Cannot open file $0: $!\n";
	while(<MYSELF>) {
		last if (/^\#TAGSTART/);
	}
	my $available = "";
	while(<MYSELF>) {
		last if (/^\#TAGEND/);
		if (/(^\d+)\)/) {
			$available = $available.$1." ";
		} 
		print if ($needsPrinting);
	}
	close(MYSELF);

	my @temp = split(/ /,$available);
	for (my $i=0;$i<=$#temp;$i++)  {
		$_ = $temp[$i] + 100;
		$available = $available." ".$_;
	}
	return $available;
}


