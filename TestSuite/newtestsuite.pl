#! /usr/bin/perl -w

=pod
// BEGIN LICENSE BLOCK
/*
Copyright ? 2008 , UT-Battelle, LLC
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
my ($executable, $inputFile,$dataFile,$testNumber,$all);  #missing input and data implementation
my $hashTable = "hashTable.txt";
my ($generateFlag,$buildFlag,$runFlag,$cmpFlag) = (1,1,1,1);	#Enable ?,?, runSingleTest, ?

GetOptions("exe=s" => \$executable, "inp=s" => \$inputFile, "dat=s" => \$dataFile, "n=i" => \$testNumber, "all" => \$all, "g=i" => \$generateFlag, "b=i" => \$buildFlag, "r=i" => \$runFlag, "c=i" => \$cmpFlag);	#Provides command options when running the script

#my @GlobalHashTable = loadHashTable($hashTable);

& selectTests();	#Starts running the script

sub selectTests
{
	if(!defined($all) && !defined($testNumber)) {		#If no test is explicitly defined
		$testNumber = askWhatTest();	#Make the user select a test from the available ones
	}
	
	if($all) {
		runAllTests(0);		#Run all tests
	} elsif($testNumber < 0) {
		runAllTests(-$testNumber);	#Run all tests starting at $testNumber
	} else {
		testSuite($testNumber);	#Run specified test
	}
}

sub askWhatTest
{
	my $available = getAvailableTests();	#Get a numerical list of the available tests
	my @temp;
	my @testArray;
	
	#Here the user selects a test to be run
	while() {
		print "Please type the number of the test you want to run\n";
		print "Available: $available\n";
		print "Default: 0 (press ENTER) ";
		chomp(my $tn = <STDIN>);	#Input test number
		$tn = 0 if($tn eq "");
		@testArray = split(/ /,$available);	#Split string into test numbers
		@temp = grep(/$tn/,@testArray);		#Search if the wanted test is available
		
		if(!@temp) {		#If user selection is not a valid test number then produce error
			print "\nERROR: An incorrect test was selected! Try again.\n\n";
		} else {
			print "OK, I'll run test number $tn\n";
			undef @temp;
			return $tn;		#If valid selection then continue to run the test
		}
	}
}

sub getAvailableTests
{
	my $available = "";
	
	open(THIS,$0) || die "ERROR: Cannot open file $0: $!\n";	#Open this file: testsuite.pl
	while(<THIS>) {
		last if(/^\#TAGEND/);		#TAGEND delimits the available tests in this document
		if(/(^\d+)\)/) {
			$available = $available.$1." ";	
		}
		print;	#Display the available tests with their descriptions provided
	}
	close(THIS);	

	my @testsArray = split(/ /,$available);	#Split the numbers of the tests
	
	for(my $i = 0;$i <= $#testsArray; $i++) {	
		$_ = $testsArray[$i] + 100;
		$available = $available." ".$_;		#Add to the list the tests with no SU(2) symmetry
	}
	
	return $available;
}

sub validateFile	#Verifies if the file given exists
{
	my ($file) = @_;

	if(-e $file) {		#Check if the file exists
		return 1;
	} else {
		print "ERROR: The input file does not exists.\n";
		return 0;
	}
}

sub executableExists
{
	my ($file) = @_;
	my $hashKey;
	
	if(validateFile($file)) {
		$hashKey = createHashKey($file);
	
		if(foundHashKey($hashKey,$hashTable)) {
			$executable = "../src/dmrg-$hashKey";
			return 1;
		} else {
			print "The hash key for $file was not found in the hash table.\n";
			return 0;
		}
	} else {
		print "The file $file does not exists.\n";
		return 0;
	}
}

sub createHashKey
{
	my ($file) = @_;
	
	open(MYHASH,"md5sum $file|") || die "Cannot open for system command: $!\n";	#Use 'md5sum' to create 128-bit hash key
	my $hashKey = substr(<MYHASH>,0,10);	#Use only the first 10 characters as the hash key
	close(MYHASH);
	
	return $hashKey;
}

sub foundHashKey	#Verify if a certain hash key is found in the output file
{
	my ($searchKey,$outfile) = @_;
	my $key;
	my $found = 0;

	open(FILE,$outfile) || die "Cannot open file $outfile: $!\n";
	while(<FILE>) {		#Iterate through the all the hash values
		chomp;
		chomp;
		$found = 1 if($_ eq $searchKey);	#Return 'true' if the hash key was matched
	}	
	close(FILE);

	return $found;
}

#-------------------------------until here is good
sub testSuite  #this function should verify for $executable, $inputFile , $dataFile
{
	my ($tn) = @_;
	
	if(defined($executable)) {
		#do stuff with  executable
		#validate the name of executable
		runSingleTest($tn) if($runFlag);
	} elsif(defined($inputFile)) {
		#do stuff with input.inp
	} elsif(defined($dataFile)) {
		#do stuff with data.txt
	} else {
		my $searchFile = "inputs/model$tn".".spec";
		
		if(!executableExists($searchFile)) {	#verify if the executable for test number already exists
			$executable = createExecutable($tn);	#creates executable
			print "Executable was created.\n";
		} else {
			print "Executable existed.\n";
		}
		
		runSingleTest($tn) if($runFlag);
	}
	
	
	#if(inHashTable($tn)) {
	#	print "Test number $tn has already been built and compiled.\n";
	#}
		
	if($cmpFlag) {
		#validateEnergy($tn);
		#validateProfile($tn);
		#S,N,C	
	}
}

sub createExecutable
{
	my ($tn) = @_;
	$tn -= 100 if($tn >= 100);
	print "Trying to create executable for test $tn...\n";
	chdir("../src");
	my $specFile = "../TestSuite/inputs/model$tn.spec";
	my $rCode = 0;
	
	if($generateFlag) {
		$rCode = system("perl configure.pl < $specFile >& /dev/null");
		die "Error creating executable (configure.pl failed)\n" if($rCode != 0);
	}
	
	if($buildFlag) {
		$rCode = system("make -f Makefile");
		die "Error creating executable (make failed)\n" if($rCode != 0);
	}
	my $hashKey = createHashKey($specFile);
	system("mv dmrg dmrg-$hashKey");
	chdir("../TestSuite");
	print "Executable was succesfully created.\n";
	
	return "../src/dmrg";	
}

sub runSingleTest
{
	my ($tn) = @_;
	print "Please wait while the test is run...\n";
	unlink("data$tn.txt");
	my $rCode = system("$executable inputs/input$tn.inp >& /dev/null ");
	die "$0: Test cannot be run to completion\n" if ($rCode != 0);
	print "The run has been completed.\n";
}

#----------------------------------until here under construction


sub runAllTests
{
	print "Run all tests starting from test #@_\n";
	my ($start) = @_;
	my $available = getAvailableTests(0);
	my @temp = split(/ /,$available);
	
	print "All tests done succesfully.\n";
}