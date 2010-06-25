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
=cut

use strict;
use Getopt::Long;

my ($executable,$testNumber,$all,$update); 
my $hashTable = "hashTable.txt";
my ($generateFlag,$buildFlag,$runFlag,$cmpFlag,$observeFlag) = (1,1,1,1,0);	#Enable driver,enable makefile, enable test run,enable validations, enable cooking of observables

GetOptions("exe=s" => \$executable, "n=i" => \$testNumber, "all" => \$all, "update" => \$update, "g=i" => \$generateFlag, "b=i" => \$buildFlag, "r=i" => \$runFlag, "c=i" => \$cmpFlag, "o=i" => \$observeFlag);	#Provides command options when running the script

if($update) {
	updateHashTable();
} 

selectTests();	#Starts running the script

#~ sub updateHashTable  #incompleete
#~ {
	#~ chdir("../src");
	#~ my @searchFiles = glob "dmrg-*";	#Put into an array all the executables found
	#~ print "Initiating update process...\n";
	#~ chdir("../TestSuite");
	#~ my (@searchKeys, $tempK, $file);
	
	#~ if(!@searchFiles) {
		#~ system("echo -n > $hashTable >& /dev/null");
		#~ die "All keys were deleted from $hashTable.\n";
	#~ }
		
	#~ foreach $file(@searchFiles) {
		#~ push(@searchKeys,substr($file,5,10));
	#~ }
		
	#~ foreach $tempK(@searchKeys) {
		#~ if(!findKey($tempK,$hashTable)) {
			#~ addKey($tempK);
		#~ }		
	#~ }
	
	#~ my @contents;
	#~ open(MYTABLE, $hashTable) || die "Cannot open file $hashTable: $!\n";
	#~ while($tempK = <MYTABLE>) {		#Iterate through the all the hash values
		#~ chomp($tempK);
		#~ chomp($tempK);
		#~ push(@contents, $tempK);
		#~ print @contents."--".@searchKeys;
		#~ #if (grep(!/$tempK/, @searchKeys)) {
			#~ #@contents = grep(!/$tempK/, @contents);
			#~ #print @contents."hola\n";
		#~ }
	#~ close(MYTABLE);
#~ }

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
	my $available = getAvailableTests(1);	#Get a numerical list of the available tests
	my @temp;
	my @testArray;

	#Here the user selects a test to be run
	while() {
		print "Please type the number of the test you want to run\n";
		print "Available: $available\n";
		print "Default: 0 (press ENTER) ";
		chomp(my $tn = <STDIN>);	#Input test number
		
		if($tn eq "") {
			$tn = 0;
		}
		
		my $searchNum = abs($tn);
		@testArray = split(/ /,$available);	#Split string into test numbers
		@temp = grep(/$searchNum/,@testArray);		#Search if the wanted test is available
		
		if(!@temp) {		#If user selection is not a valid test number then produce error
			print "\nERROR: An incorrect test was selected! Please try again.\n\n";
		} else {
			print "Preparing to run test number $tn\n";
			undef @temp;
			return $tn;		#If valid selection then continue to run the test
		}
	}
}

sub getAvailableTests
{
	my ($pFLAG) = @_;
	my $available = "";
	my $descriptionFile = "inputs/descriptions.txt";
	
	open(FILE,$descriptionFile) || die "ERROR: Cannot open file $descriptionFile: $!\n";	#Open this file: testsuite.pl
	while(<FILE>) {
		last if(/^\#TAGEND/);		#TAGEND delimits the available tests in this document
		if(/(^\d+)\)/) {
			$available = $available.$1." ";	
		}
		print if($pFLAG);	#Display the available tests with their descriptions provided
	}
	close(FILE);	

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
		print "ERROR: The file $file does not exists.\n";
		return 0;
	}
}

sub removeKey
{
	my ($mkey,$table) = @_;
	my @contents;
	my $line;
	
	open FILE, "<$table";
	while($line = <FILE>) {
		last if($line eq "" || $line eq "\n");
		next if(substr($line,0,10) eq $mkey) ;
		push(@contents,$line);
	}
	close FILE;
	shift(@contents);
	open (OUT, ">$table");
	foreach $line(@contents) {
		print  OUT $line;
	}
	print OUT "\n";
	close (OUT);
	
	print "The key $mkey has been removed.\n";
}

sub executableExists
{
	my ($tn, $model) = @_;
	my $modelKey = createHashKey($model);

	if(findKey($modelKey,$hashTable)) {
		if(!findTn($modelKey,$tn,$hashTable)) {
			addTn($modelKey,$tn,$hashTable);
		}
		
		$executable = "../src/dmrg-$modelKey";
		return 1;
		
	} else {
		print "The hash key for $model was not found in the hash table.\n";
		return 0;
	}
}

sub addTn
{
	my ($key,$tn,$table) = @_;
	my @contents;
	my $line;
	my $old;
	
	open (FILE, "<$table");
	while($line = <FILE>) {
		last if($line eq "" || $line eq "\n");
		if(substr($line,0,10) eq $key) {
			chomp($line);
			$old = $line;
		} else {
			push(@contents,$line);
		}
	}
	close (FILE);
	
	my $new = $old." ".$tn;
	push(@contents,$new);
	
	open (OUT, ">$table");
	foreach $line(@contents) {
		print  OUT $line;
	}
	print OUT "\n";
	close (OUT);
	
	print "Test number was added: $new.\n";
}

sub addKey
{
	my ($mkey,$tn,$table) = @_;
	
	open FILE, ">>$table";
	print FILE "$mkey - $tn\n";
	close FILE;
	
	print "New key: $mkey - $tn\n";	
}

sub findKey
{
	my ($mkey,$table) = @_;
	open FILE, "<$table" or die "Can't open $table: $!";
	while (my $line = <FILE>) {
		if (substr($line,0,10) eq $mkey) {
			close FILE;
#			print "Model key found.\n";
			return 1;
		}
	}
	close FILE;
	
#	print "Model key not found.\n";
	return 0;	
}


sub findTn
{
	my ($key,$tn,$table) = @_;
	my $line;
	my $found = 0;
	
	open FILE, $table;
	do { $line = <FILE> } until substr($line,0,10) eq $key || eof;
	close FILE;

	my @temp = split(" ",$line);
	shift(@temp);
	shift(@temp);
	
	if(!grep {$_ == $tn} @temp) {
		print "Test number $tn has not been compiled.\n";
	} else {
		$found = 1;
		print "Test number $tn is already compiled.\n";
	}
	
	return $found;
}

sub createHashKey
{
	my ($file) = @_;
	
	open(MYHASH,"md5sum $file|") || die "Cannot open for system command: $!\n";	#Use 'md5sum' to create 128-bit hash key
	my $hashKey = substr(<MYHASH>,0,10);	#Use only the first 10 characters as the hash key
	close(MYHASH);
	
	return $hashKey;
}

sub testSuite  #this function should verify for $executable, $inputFile , $dataFile
{
	my ($tn) = @_;
	my $inputFile = "inputs/input$tn.inp";
	
	my $temp = $tn;
	$temp -= 100 if($tn >= 100);
	my $modelFile = "inputs/model$temp".".spec";
	die "An input file is missing: $inputFile\n" if(!validateFile($inputFile));
	
	if(defined($executable)) {
		die "Missing executable file: $executable" if(!validateFile($executable));
		
		#think about the hash and hashTable
	
	} else {
		if(validateFile($modelFile)) {
			if(!executableExists($tn,$modelFile)) {	#verify if the executable for test number already exists
				$executable = createExecutable($tn);	#creates executable
			} else {
				print "Retrieving existing executable.\n";
			}
		} else {   #die
		die "Execution aborted: $modelFile is missing.\n";
		}
	}
	
	runSingleTest($tn,$executable) if($runFlag);
	undef $executable;	
	
	if($cmpFlag) {
		#validateEnergy($tn);
		#validateProfile($tn);
		#S,N,C	
	}
}

sub createExecutable
{
	my ($tn) = @_;
	chdir("../src");

	$tn -= 100 if($tn >= 100);
	print "Creating executable for test $tn...\n";
	print "***********************************\n";
	
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
	print "***********************************\n";
	my $specKey = createHashKey($specFile);
	$rCode = system("mv dmrg dmrg-$specKey");
	die "Error renaming the executable file\n" if($rCode != 0);
	chdir("../TestSuite");
	addKey($specKey,$tn,$hashTable);
	print "Executable was succesfully created.\n";
	
	return "../src/dmrg-$specKey";	
}

sub runSingleTest
{
	my ($tn,$executable) = @_;
	print "Please wait while the test is run...\n";
	unlink("data$tn.txt");
	my $rCode = system("$executable inputs/input$tn.inp >& /dev/null ");
	die "$0: Test cannot be run to completion\n" if ($rCode != 0);
	print "The run has been completed.\n";
}

sub runAllTests
{
	my ($start) = @_;
	print "Running all tests starting from test #$start\n";
	
	my $available = getAvailableTests(0);
	my @temp = split(/ /,$available);
	
	for (my $i=0;$i<$#temp+1;$i++) {
		next if ($temp[$i] eq "");
		next if ($temp[$i]<$start);
		testSuite($temp[$i]);
	}
	print "All tests done succesfully.\n";
}