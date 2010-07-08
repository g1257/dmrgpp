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

my ($executable,$observable,$testNumber,$all,$update,$verbose,$noModel) = (undef,undef,undef,0,0,0,0); 
my ($cmpFlag,$observeFlag,$rmFlag) = (1,0,0);

GetOptions("exe=s" => \$executable, "obs=s" => \$observable, "n=i" => \$testNumber, "all" => \$all, "u" => \$update,
"c" => \$cmpFlag, "o" => \$observeFlag, "rm" => \$rmFlag, "v" => \$verbose, "m" => \$noModel);	#Provides command options when running the script

my $srcDir = "../src/";
my $testDir = "../TestSuite/";
my $oraclesDir = $testDir."oracles/";		#Specifies the output folder for the oracles results
my $resultsDir = $testDir."results/";
my $inputsDir = $testDir."inputs/";
my $hashTable = $testDir."hashTable.txt";


runScript();

sub updateHashTable  #incompleete
{
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
}

sub runScript
{
	if($update) {
		updateHashTable();
	} 

	if(validateFile($hashTable) && validateDirectory($srcDir) && validateDirectory($testDir) && validateDirectory($oraclesDir) && validateDirectory($resultsDir) && validateDirectory($inputsDir)) {
		if(!($all) && !defined($testNumber)) {		#If no test is explicitly defined
			$testNumber = selectTest();	#Make the user select a test from the available ones
		}
		
		if($all) {
			runAllTests(0);		#Run all tests
		} elsif($testNumber < 0) {
			runAllTests(-$testNumber);	#Run all tests starting at $testNumber
		} else {
			runTest($testNumber);	#Run specified test
		}
	}
}

sub validateDirectory	
{
	my ($dir) = @_;
	
	if(!-d $dir) {
		if($dir eq $oraclesDir || $dir eq $resultsDir) {
			my $rCode = system("mkdir $dir");
			die "Error making $dir: $!\n" if($rCode);
			print "Directory created: $dir\n";
		} else {
			print "Error searching for $dir: $!\n";
			return 0;
		}
	}
	
	return 1;
}

sub selectTest
{
	my $available = getAvailableTests();	#Get a numerical list of the available tests
	my @temp;
	my @testArray;
	my $searchNum;
	#Here the user selects a test to be run
	while() {
		print "Type the number of the test you want to run.\n";
		print "Available tests: $available\n";
		print "Default is 0 (press ENTER): ";
		chomp(my $tn = <STDIN>);	#Input test number
		$tn = 0 if($tn eq "");
		
		$searchNum = abs($tn);
		@testArray = split(/ /,$available);	#Split string into test numbers
		@temp = grep(/$searchNum/,@testArray);		#Search if the wanted test is available
		
		if(!@temp) {		#If user selection is not a valid test number then produce error
			print "Error: An incorrect test was selected! Try again.\n";
		} else {
			return $tn;		#If valid selection then continue to run the test
		}
	}
}

sub getAvailableTests
{
	my $available = "";
	my $descriptionFile = $inputsDir."descriptions.txt";
	
	open(FILE,$descriptionFile) || die "Error opening $descriptionFile: $!\n";	#Open this file: testsuite.pl
	while(<FILE>) {
		last if(/^\#TAGEND/);		#TAGEND delimits the available tests in this document
		if(/(^\d+)\)/) {
			$available = $available.$1." ";	
		}
		print if($verbose);	#Display the available tests with their descriptions provided
	}
	close(FILE);	

	my @testsArray = split(/ /,$available);	#Split the numbers of the tests
	my $temp;
	
	for(my $i = 0;$i <= $#testsArray; $i++) {	
		$temp = $testsArray[$i] + 100;
		$available = $available." ".$temp;		#Add to the list the tests with no SU(2) symmetry
	}
	
	return $available;
}

sub validateFile	#Verifies if the file given exists
{
	my ($file) = @_;

	if(-e $file) {		#Check if the file exists
		return 1;
	} elsif($file eq $hashTable) {
		my $rCode = system("touch $hashTable");
		print "Error creating hash table: $!\n";
		return 1;
	} else {
		#print "Error: The file $file does not exists.\n";
		return 0;
	}
}

sub removeKey
{
	my ($mkey) = @_;
	my @contents;
	my $line;
	
	open FILE, "<$hashTable";
	while($line = <FILE>) {
		last if($line eq "" || $line eq "\n");
		next if(substr($line,0,10) eq $mkey) ;
		push(@contents,$line);
	}
	close FILE;
	shift(@contents);
	open (OUT, ">$hashTable");
	foreach $line(@contents) {
		print  OUT $line;
	}
	print OUT "\n";
	close (OUT);
	
	print "The key $mkey has been removed.\n";
}

sub removeFiles
{
	my ($tn) = @_;
	my $files = "Makefile* main* observe freeSystem* dmrg gmon.out observe.* input.* raw$tn.txt gmon.out";
	my $rCode;
	my @temp = split(/ /,$files);
	
	foreach my $f (@temp) {
		$rCode = system("rm $f") if(validateFile("$f"));
	}
	
	chdir("$srcDir");
	foreach my $f (@temp) {
		$rCode = system("rm $f") if(validateFile("$f"));
	}
	chdir("$testDir");
	
	print "All temporary files were removed.\n";
}

sub executableExists
{
	my ($key) = @_;
	
	if(findKey($key)) {
		print "Retrieving existing executable...\n";
		$executable = $srcDir."dmrg-$key";
		die "Error retrieving $executable: $!\n" if(!validateFile($executable));
		return 1;
	} else {
		print "The executable was not found.\n";
		return 0;
	}
}

sub observableExists
{
	my ($key) = @_;
	
	if(findKey($key)) {
		print "Retrieving existing observable...\n";
		$observable = $srcDir."observe-$key";
		die "Error retrieving $observable: $!\n" if(!validateFile($observable));
		return 1;
	} else {
		print "The observable was not found.\n";
		return 0;
	}
}

sub addTn
{
	my ($key,$tn) = @_;
	my @contents;
	my $line;
	my $old;
	
	open (FILE, "<$hashTable");
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
	
	open (OUT, ">$hashTable");
	foreach $line(@contents) {
		print  OUT $line;
	}
	print OUT "\n";
	close (OUT);
	
	print "Test number was added: $new.\n";
}

sub addKey
{
	my ($mkey,$tn) = @_;
	
	open FILE, ">>$hashTable";
	print FILE "$mkey - $tn\n";
	close FILE;
	
	print "New key: $mkey - $tn\n";	
}

sub findKey
{
	my ($mkey) = @_;
	open FILE, "<$hashTable" or die "Can't open $hashTable: $!";
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
	my ($key,$tn) = @_;
	my $line;
	my $found = 0;
	
	open FILE, $hashTable;
	do { $line = <FILE> } until substr($line,0,10) eq $key || eof;
	close FILE;

	my @temp = split(" ",$line);
	shift(@temp);
	shift(@temp);
	
	if(!grep {$_ == $tn} @temp) {
		print "Test number $tn was not found.\n";
	} else {
		$found = 1;
		print "Test number $tn was found.\n";
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

sub runTest  #this function should verify for $executable, $inputFile , $dataFile
{
	my ($tn) = @_;
	my $temp = $tn;
	$temp -= 100 if($tn >= 100);
	my $specFile = $testDir.$inputsDir."model$temp.spec";
	my $rCode;
	my $addFlag = 0;
	my $configFile = $srcDir."configure.pl";
	
	die "Error: Model file is missing.\n" if(!validateFile($specFile));
	my $specKey = createHashKey($specFile);
	
	if(defined($executable)){
		die "Missing executable file: $executable" if(!validateFile($executable));
		die "Missing observable file: $observable" if(!validateFile($observable));
	}elsif($noModel) {
		chdir("$srcDir");
		$rCode = system("./$configFile");
		die "Error creating executable $configFile: $!\n" if($rCode);
		$rCode = system("make -f Makefile");
		die "Error with the Makefile: $!\n" if($rCode);
		$executable = $srcDir."dmrg";
		chdir("$testDir");
	}elsif(!executableExists($specKey)) {
		$executable = createExecutable($tn,$specFile,$specKey,$configFile);
		$addFlag = 1;
	} else {
		if(!findTn($specKey,$tn) ){
			addTn($specKey,$tn);
		}
	}
	
	if($observeFlag) {
		if(!observableExists($specKey)) {
			$observable = createObservable($tn,$specKey);
		}
	}
	addKey($specKey,$tn) if(!$noModel && $addFlag);
	print "Please wait while the test is run...\n";
	
	chdir("$srcDir");
	my $inputFile = $testDir.$inputsDir."input$tn.inp";	
	$rCode = system("$executable $inputFile >& /dev/null ");
	die "$0: Test cannot be run to completion\n" if ($rCode);
	print "The run has been completed.\n";
	my $profFile = $testDir.$resultsDir."prof$tn.txt";
	profile($profFile);
	chdir("$testDir");
	
	my $dataFile = $srcDir."data$tn.txt";
	my $rawFile = $testDir."raw$tn.txt";
	my $energyOut = $testDir.$resultsDir."e$tn.txt";
	my $tstFile = $srcDir."tst$tn.txt";
	my $envStack = $srcDir."EnvironStackdata$tn.txt";
	my $sysStack = $srcDir."SystemStackdata$tn.txt";
	
	observables($tn,$inputFile,$dataFile,$rawFile) if($observeFlag);
	extractEnergy($dataFile,$energyOut);
	moveFiles($dataFile,$tstFile,$envStack,$sysStack);
	removeFiles($tn) if($rmFlag);
	print "Completed running test $tn.\n";

	if($cmpFlag) {
		validateEnergy($tn);
		validateProfile($tn);
		validateSNC($tn);	
	}
	
	print "******END OF TEST ".$tn."********\n";
	
	undef $executable;
	undef $observable;
	undef $tn;
	
}

sub moveFiles
{
	my (@files) = @_;
	my $rCode;
	
	foreach my $f (@files) {
		$rCode = system("mv $f ".$testDir.$oraclesDir) if(validateFile("$f"));
	}
	
	chdir("$srcDir");
	foreach my $f (@files) {
		$rCode = system("mv $f ".$testDir.$oraclesDir) if(validateFile("$f"));
	}
	chdir("$testDir");
	
	print "All results were succesfully stored.\n";
}

sub validateEnergy
{
	
	
}

sub validateProfile
{
	
	
}

sub validateSNC
{
	
	
}

sub createObservable
{
	my ($tn,$specKey) = @_;
	my $tempExe;
	
	print "Creating observable for test $tn...\n";
	chdir("$srcDir");
	my $rCode;
	
	if($noModel) {
		$rCode = system("make observe -f Makefile");
		die "Error with the Makefile: $!\n" if($rCode);
		$tempExe = "observe";
	} else {
		$rCode = system("make observe -f Makefile");
		die "Error with the Makefile: $!\n" if($rCode);
		$tempExe = "observe-$specKey";
		$rCode = system("mv observe $tempExe");
		die "Error renaming the observable file\n" if($rCode);
	}
	chdir("$testDir");

	print "Observable was succesfully created.\n";
	
	return $srcDir.$tempExe;
}

sub createExecutable
{
	my ($tn,$spec,$specKey,$configFile) = @_;
	my $tempExe;
	
	print "Creating executable for test $tn...\n";
	chdir("$srcDir");
	my $rCode;
	
	$rCode = system("./$configFile < $spec >& /dev/null");
	die "Error creating executable $configFile: $!\n" if($rCode);
	$rCode = system("make -f Makefile");
	die "Error with the Makefile: $!\n" if($rCode);
	$tempExe = "dmrg-$specKey";
	$rCode = system("mv dmrg $tempExe");
	die "Error renaming the executable file\n" if($rCode);
	chdir("$testDir");

	print "Executable was succesfully created.\n";
	
	return $srcDir.$tempExe;	
}

sub observables
{
	my ($tn,$input,$data,$raw) = @_;
	my $rCode;
	
	$rCode = system($srcDir."$observable $input $data > $raw"); 
	die "Error running observables: $!\n" if($rCode);

	my $cOut = $testDir.$resultsDir."operatorC$tn.txt";
	my $nOut = $testDir.$resultsDir."operatorN$tn.txt";
	my $sOut = $testDir.$resultsDir."operatorS$tn.txt";
	
	extractOperatorC($raw,$cOut);
	extractOperatorN($raw,$nOut);
	extractOperatorS($raw,$sOut);
}


sub extractEnergy
{
	my ($data,$energyOut) = @_;
        my $rCode = system("grep Energy $data > $energyOut");
	die "Error extracting energy values: $!\n" if($rCode);
        print "Energy extraction was succesful!\n";
}

sub profile
{
	my ($prof) = @_;
	my $rCode = system("gprof $executable > $prof"); 
	die "Error with gprof: $!\n" if($rCode); 
	 print "Profiling was succesful!\n";
}

sub extractOperatorC
{
	my ($raw,$cOut) = @_;
	my $line;
	my $opC;
	open(INFILE,"<$raw") || die "Error opening file: $!\n";
		while($line = <INFILE>) {
			if($line =~ /^OperatorC/) {
				$opC = $line;
				$line = <INFILE>;
				$opC = $opC.$line;
				my @temp1 = split(/ /,$line);

				for(my $i = 0; $i < $temp1[0]; $i++) {
					$line = <INFILE>;
					$opC = $opC.$line;
				}
			}
		}
	close(INFILE);
	
	open (OUTFILE, ">$cOut");
	print OUTFILE $opC;
	close (OUTFILE);
	print "OperatorC extraction was succesful!\n";
}


sub extractOperatorN
{
	my ($raw,$nOut) = @_;
	my $line;
	my $opN;
	open(INFILE,"<$raw") || die "Error opening file: $!\n";	#Open this file: testsuite.pl
		while($line = <INFILE>) {
			
			if($line =~ /^OperatorN/) {
				$opN = $line;
				$line = <INFILE>;
				$opN = $opN.$line;
				my @temp1 = split(/ /,$line);

				for(my $i = 0; $i < $temp1[0]; $i++) {
					$line = <INFILE>;
					$opN = $opN.$line;
				}
				
			}
		}
	close(INFILE);
	
	open (OUTFILE, ">$nOut");
	print OUTFILE $opN;
	close (OUTFILE);
	print "OperatorN extraction was succesful!\n";
}


sub extractOperatorS
{
	my ($raw,$sOut) = @_;
	my $line;
	my $opS;
	open(INFILE,"<$raw") || die "Error: Cannot open file: $!\n";
		while($line = <INFILE>) {
			
			if($line =~ /^OperatorS/) {
				$opS = $line;
				$line = <INFILE>;
				$opS = $opS.$line;
				my @temp1 = split(/ /,$line);

				for(my $i = 0; $i < $temp1[0]; $i++) {
					$line = <INFILE>;
					$opS = $opS.$line;
				}
				
			}
		}
	close(INFILE);
	
	open (OUTFILE, ">$sOut");
	print OUTFILE $opS;
	print "OperatorSz extraction was succesful!\n";
}

sub runAllTests
{
	my ($start) = @_;
	print "Preparing to run all tests starting from test $start..\n";
	
	my $available = getAvailableTests(0);
	my @temp = split(/ /,$available);
	my @notFunctionalTests = (4,24,60,104,105,106,107,108,124,125,141,160);
	
	for (my $i=0;$i<=$#temp;$i++) {
		next if ($temp[$i] eq "");
		next if ($temp[$i]<$start);
		next if(grep {$_ eq $temp[$i]}@notFunctionalTests);
		runTest($temp[$i]);
	}
	print "All tests done succesfully.\n";
}