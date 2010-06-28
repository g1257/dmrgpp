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

my ($executable,$testNumber,$all,$update,$verbose,$noModel) = (undef,undef,0,0,0,0); 
my ($generateFlag,$buildFlag,$runFlag,$cmpFlag,$observeFlag,$rmFlag) = (1,1,1,1,0,0);	#Enable driver,enable makefile, enable test run,enable validations, enable cooking of observables

GetOptions("x=s" => \$executable, "n=i" => \$testNumber, "all" => \$all, "u" => \$update,
"g" => \$generateFlag, "b" => \$buildFlag, "r" => \$runFlag, "c" => \$cmpFlag, "o" => \$observeFlag,
"rm" => \$rmFlag, "v" => \$verbose, "m" => \$noModel);	#Provides command options when running the script

my $oraclesDir = "oracles/";		#Specifies the output folder for the oracles results
my $srcDir = "../src/";
my $testDir = "../TestSuite/";
my $hashTable = "hashTable.txt";
my $resultsDir = "results/";


selectTest();

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

sub selectTest
{
	if($update) {
		updateHashTable();
	} 

	if(validateFile($hashTable) && validateDirectory($srcDir) && validateDirectory($testDir) && validateDirectory($oraclesDir) && validateDirectory($resultsDir)) {
		if(!($all) && !defined($testNumber)) {		#If no test is explicitly defined
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
			die "Error searching for $dir: $!\n";
		}
	}
	
	return 1;
}

sub askWhatTest
{
	my $available = getAvailableTests();	#Get a numerical list of the available tests
	my @temp;
	my @testArray;
	my $searchNum;
	#Here the user selects a test to be run
	while() {
		print "Please type the number of the test you want to run\n";
		print "Available: $available\n";
		print "Default: 0 (press ENTER) ";
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
	my $descriptionFile = "inputs/descriptions.txt";
	
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
		print "Error: The file $file does not exists.\n";
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

sub removeFiles
{
	my ($tn) = @_;
	my $rmTest = "raw$tn.txt gmon.out";
	my $rmSrc = "Makefile* main* freeSystem* dmrg gmon.out observe.* input.*";
	my $rCode;
	my $tst = 0;
	
	if($observeFlag) {
		$rCode = system("rm $rmTest");
		#die "Error removing file: $!\n" if($rCode);
		$rmSrc = $rmSrc." observe";
	}
	
	chdir("$srcDir");
	$rCode = system("rm $rmSrc");
	#die "Error removing file: $!\n" if($rCode);
	chdir("$testDir");
	print "All temporary files were removed.\n";
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
	my $inputFile = $testDir."inputs/input$tn.inp";
	my $temp = $tn;
	$temp -= 100 if($tn >= 100);
	my $specFile = $testDir."inputs/model$temp.spec";
	
	exit() if(!validateFile($inputFile) || !validateFile($specFile));
	print "Preparing to run test $tn...\n";
	
	if(defined($executable)) {
		die "Missing executable file: $executable" if(!validateFile($executable));
		
		#think about the hash and hashTable
	
	} elsif(!executableExists($tn,$specFile)) {	#verify if the executable for test number already exists
				$executable = createExecutable($tn,$specFile);	#creates executable
	} else {
		print "Retrieving existing executable.\n";
	}
	
	runSingleTest($tn,$executable,$specFile,$inputFile) if($runFlag);
	#--------------------
	system("mv data$tn.txt results/");
	if($cmpFlag) {
		validateEnergy($tn);
		validateProfile($tn);
		validateSNC($tn);	
	}
	
	undef $executable;
	undef $tn;
	print "******END OF TEST $tn********\n";
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


sub createExecutable
{
	my ($tn,$spec) = @_;
	
	chdir("$srcDir");
	$tn -= 100 if($tn >= 100);
	print "Creating executable for test $tn...\n";
	my $rCode;
	
	print "***********************************\n";
	if($generateFlag) {
		$rCode = system("./configure.pl < $spec >& /dev/null");
		die "Error creating executable (configure.pl failed)\n" if($rCode);
	}

	if($buildFlag) {
		$rCode = system("make -f Makefile");
		die "Error creating executable (make failed)\n" if($rCode);
	}
	print "***********************************\n";
	my $specKey = createHashKey($spec);
	$rCode = system("mv dmrg dmrg-$specKey");
	die "Error renaming the executable file\n" if($rCode);
	chdir("$testDir");
	
	addKey($specKey,$tn,$hashTable);
	print "Executable was succesfully created.\n";
	
	return $srcDir."dmrg-$specKey";	
}

sub runSingleTest
{
	my ($tn,$executable,$spec,$input) = @_;
	print "Please wait while the test is run...\n";
	my $profFile = $testDir.$resultsDir."prof$tn.txt";
	my $rCode = system("$executable $input >& /dev/null ");
	die "$0: Test cannot be run to completion\n" if ($rCode);
	print "The run has been completed.\n";
	profile($executable,$profFile);
	my $hashKey = createHashKey($spec);
	my $dataFile = $srcDir."data$tn.txt";
	my $rawFile = "raw$tn.txt";
	my $energyOut = $resultsDir."e$tn.txt";
	
	
	extractEnergy($dataFile,$energyOut);
	observables($tn,$input,$dataFile,$rawFile) if($observeFlag);
	$rCode = system("mv $dataFile $resultsDir");
	die "Error moving $dataFile: $!\n" if($rCode);
	removeFiles($tn) if($rmFlag);
	print "Completed running test $tn.\n";
}

sub observables
{
	my ($tn,$input,$data,$raw) = @_;

	chdir("$srcDir");
	my $rCode = system("make observe -f Makefile");
	die "Error with making observe: $!\n" if($rCode);
	chdir("$testDir");
	
	$rCode = system($srcDir."observe $input $data > $raw"); 
	die "Error running observables: $!\n" if($rCode);

	my $cOut = $resultsDir."operatorC$tn.txt";
	my $nOut = $resultsDir."operatorN$tn.txt";
	my $sOut = $resultsDir."operatorS$tn.txt";
	
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
	my ($exe,$prof) = @_;
	my $rCode = system("gprof $exe > $prof"); 
	die "Error with gprof: $!\n" if($rCode); 
	 print "Profiling was succesful!\n";
}

#~ sub cookRaw
#~ {
	#~ my ($tn,$input,$data,$output,$hashKey) = @_;
	
	#~ chdir("../src/");
	#~ my $rCode = system("make observe -f Makefile");
	#~ die "Error make observe.\n" if($rCode);
	#~ $rCode = system("mv observe observe-$hashKey");
	#~ die "Error renaming the executable file\n" if($rCode != 0);
	#~ chdir("../TestSuite/");
	
	#~ $rCode = system("../src/observe-$hashKey $input $data > $output"); 
	#~ die "Error execute observe.\n" if($rCode);
	
	#~ my $line;
	#~ my $opC;
	#~ my $opN;
	#~ my $opS;

	#~ open(FILE,"<$output") || die "ERROR: Cannot open file: $!\n";	#Open this file: testsuite.pl
		#~ while($line = <FILE>) {
			
			#~ if($line =~ /^OperatorC/) {
				#~ $opC = $line;
				#~ $line = <FILE>;
				#~ $opC = $opC.$line;
				#~ my @temp1 = split(/ /,$line);

				#~ for(my $i = 0; $i < $temp1[0]; $i++) {
					#~ $line = <FILE>;
					#~ $opC = $opC.$line;
				#~ }
				
			#~ }
			
			#~ if($line =~ /^OperatorN/) {
				#~ $opN = $line;
				#~ $line = <FILE>;
				#~ $opN = $opN.$line;
				#~ my @temp2 = split(/ /,$line);
				#~ for(my $i = 0; $i < $temp2[0]; $i++) {
					#~ $line = <FILE>;
					#~ $opN = $opN.$line;
				#~ }
				
			#~ }
		
			#~ if($line =~ /^OperatorS/) {
				#~ $opS = $line;
				#~ $line = <FILE>;
				#~ $opS = $opS.$line;
				#~ my @temp3 = split(/ /,$line);
				#~ for(my $i = 0; $i < $temp3[0]; $i++) {
					#~ $line = <FILE>;
					#~ $opS = $opS.$line;
				#~ }
				
			#~ }
			
	#~ }
	#~ close(FILE);	
	
	#~ open (FILE, ">results/operatorC$tn.txt");
	#~ print FILE $opC;
	#~ close (FILE);
	
	#~ open (FILE, ">results/operatorN$tn.txt");
	#~ print FILE $opN;
	#~ close (FILE);
	
	#~ open (FILE, ">results/operatorS$tn.txt");
	#~ print FILE $opS;
	#~ close (FILE);	
#~ }

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
	close (OUTFILE);	
	print "OperatorSz extraction was succesful!\n";
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