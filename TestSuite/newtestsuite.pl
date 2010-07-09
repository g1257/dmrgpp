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
use Cwd 'abs_path';
use File::Basename;

my ($executable,$testNum,$all) = (undef,undef,undef); 
my ($rmFlag,$update,$verbose,$noModel,$force) = (0,0,0,0,0);

GetOptions("x=s" => \$executable, "n=i" => \$testNum, "a" => \$all, "u" => \$update, 
"r" => \$rmFlag, "v" => \$verbose, "m" => \$noModel, "f" => \$force);

my ($testDir, $srcDir);
my $PATH = $testDir = $srcDir = abs_path($0);
my $filename = basename($PATH);		

$testDir =~ s/$filename$//;
$srcDir =~ s/TestSuite.*/src\//;
my $oraclesDir = $testDir."oracles/";	
my $resultsDir = $testDir."results/";
my $inputsDir = $testDir."inputs/";

my $hashTable = $testDir."hashTable.txt";
my %hash;
my $specKey;

&startProgram();

sub startProgram{
	if($update) {
		updateHashTable();
	} elsif(validateDirectory($srcDir) && validateDirectory($testDir) && validateFile($hashTable)) {
		loadHash() if(!$noModel);

		if(!($all) && !defined($testNum)) {		#If no test is explicitly defined
			selectTest();	#Make the user select a test from the available ones
		}elsif($all) {
			runAllTests(0);		#Run all tests
		}
		
		if($testNum < 0) {
			runAllTests(-$testNum) if(validateTest());	#Run all tests starting at $testNumber
		} else {
			testSuite() if(validateTest());	#Run specified test
		}
	} else {
		die "Missing source and/or test suite directories.\n";
	}
}

#~ sub updateHashTable
#~ {

#~ }

sub validateDirectory	
{
	my ($dir) = @_;
	
	if(-d $dir) {		
		return 1;
	} elsif($dir eq $resultsDir){
		system("mkdir $dir");
		print "Directory created: $dir\n";
	}
	
	return 0;
}

sub validateTest
{
	my $searchNum = abs($testNum);
	my @found = grep(/$searchNum/, split(/ /,getAvailableTests()));	
	
	if(!@found) {	
		print "\nError: An incorrect test was selected! Try again.\n\n";
		return 0;
	}

	return 1;
}

sub selectTest
{
	while() {
		print "Type the number of the test you want to run.\n";
		print "Available tests: ".getAvailableTests()."\n";
		print "Default is 0 (press ENTER): ";
		chomp($testNum = <STDIN>);	
		$testNum = 0 if($testNum eq "");
		last if(validateTest());
	}
}

sub getAvailableTests
{
	my $available = "";
	my $descriptionFile = $inputsDir."descriptions.txt";
	
	open(FILE,$descriptionFile) || die "Error opening $descriptionFile: $!\n";	
	while(<FILE>) {
		last if(/^\#TAGEND/);		
		if(/(^\d+)\)/) {
			$available .= "$1 ";	
		}
		print if($verbose);	
	}
	close(FILE);	

	my @testsArray = split(/ /,$available);	
	my $temp;
	
	for(my $i = 0;$i <= $#testsArray; $i++) {	
		$temp = $testsArray[$i] + 100;
		$available .= " $temp";		
	}
	
	return $available;
}

sub validateFile	
{
	my ($file) = @_;

	if(-e $file) {		
		return 1;
	} elsif($file eq $hashTable) {
		system("touch $hashTable >& /dev/null");
		print "File created: $hashTable\n";
		return 0;
	} else {
		return 0;
	}
} 

sub loadHash
{
	open FILE, "<$hashTable" || die "Could not open file: $hashTable\n";
	while(<FILE>) {
		chomp;
		last if($_ eq "");
		my ($key, $values) = split(/ : /);
		my @values = split(/, /, $values);
		$hash{$key} = [@values];
	}
	close FILE;
}

sub findKey
{
	my ($key) = @_;
	
	foreach my $k (keys %hash) {
		return 1 if($key eq $k);
	}
	
	return 0;
}

sub findTn
{
	my ($key) = @_;
	
	foreach my $k (keys %hash) {
		if($key eq $k) {
			if($testNum ~~ @{$hash{$key}}) {
				print "Hash table is up-to-date.\n";
				return 1;
			} 
		}
	}

	return 0;
}

sub addTn
{
	my ($key) = @_;
	
	print "Updating hash table...\n";
	
	foreach my $k (keys %hash) {
		if($key eq $k) {
			push @{$hash{$k}}, $testNum;
		}
	}
}

#~ sub removeTns
#~ {
	#~ my (@delTn) = @_;
	
	#~ print "Test numbers removed:\n";	
	#~ while(@delTn) {
		#~ my $key = shift(@delTn);
		#~ my $tn = shift(@delTn);
		#~ my @array = @{$hash{$key}};
		#~ @array = grep { !$tn } @array;
		#~ $hash{$key} = \@array;
		#~ print "$tn\n";
	#~ }
#~ }

sub addKey
{
	my ($key) = @_;
	
	print "Updating hash table...\n";
	$hash{$key} = [$testNum];
}

#~ sub removeKeys
#~ {
	#~ my (@delKeys) = @_;
	
	#~ print "Keys removed:\n";
	#~ foreach my $key (@delKeys) {
		#~ delete $hash{$key};
		#~ print "$key\n";
	#~ }
#~ }

sub saveHash
{
	open FILE, ">$hashTable" || die "Could not open file: $hashTable\n";
	foreach my $key (sort keys %hash) {
		print FILE "$key : ".join(', ', sort @{$hash{$key}})."\n";
	}
	close FILE;
}

sub createManualExecutable
{
	my ($config) = @_;
	my $arg1 = "./$config";
	my $arg2 = "make -f Makefile >& /dev/null";
	grep {s/>.*//} $arg2 if($verbose);
	
	print "Creating executable for current test...\n";
	
	chdir($srcDir);
	system($arg1);
	system($arg2); 
	chdir($testDir);
	
	print "Executable was succesfully created.\n";
	
	return $srcDir."dmrg";
}

sub createAutoExecutable
{
	my ($spec,$specKey,$configFile) = @_;
	my $arg1 = "./$configFile < $spec >& /dev/null";
	my $arg2 = "make -f Makefile >& /dev/null";
	my $exe = "dmrg-$specKey";
	
	grep {s/>.*//} $arg1 if($verbose);
	grep {s/>.*//} $arg2 if($verbose);
	
	print "Creating executable for Test $testNum...\n";
	
	chdir($srcDir);
	system($arg1);
	system($arg2);
	rename("dmrg", $exe);
	chdir($testDir);
	
	print "Executable was succesfully created.\n";
	
	return $srcDir.$exe;	
}

sub runTest
{
	my ($input) = @_;
	my $arg = "$executable $input >& /dev/null";
	
	grep {s/&//} $arg if($verbose);
	
	print "Please wait while the test is run...\n";
	chdir($srcDir);
	system($arg);
	chdir($testDir);
	print "The run has been completed.\n";	
}

sub testSuite
{
	my $temp = $testNum;
	$temp -= 100 if($testNum >= 100);
	my $specFile = $inputsDir."model$temp.spec";
	my $configFile = "configure.pl";
	my $inputFile = $inputsDir."input$testNum.inp";
	my $postFile = $inputsDir."postProcessing$testNum.txt";
	my $postProcLib = $inputsDir."postProcessingLibrary.txt";
	
	$specKey = substr(`md5sum $specFile`,0,10);
	
	if(defined($executable)){
		print "Using executable supplied.\n";
		$executable = $srcDir.$executable;
		die "Missing executable file: $executable" if(!validateFile($executable));
	}elsif($noModel) {
		$executable = createManualExecutable($configFile);
	}elsif(!findKey($specKey) || $force) {
		$executable = createAutoExecutable($specFile,$specKey,$configFile);
		addKey($specKey);
	} else {
		print "Retrieving existing executable...\n";
		$executable = $srcDir."dmrg-$specKey";
		validateFile($executable) || die "Error retrieving $executable: $!\n";
	
		if(!findTn($specKey) ){
			addTn($specKey);
		}
	}
	
	saveHash();
	#runTest($inputFile);
	if(validateFile($postFile)) {
		my @analyses = extractAnalyses($postFile) ;
		(@analyses) ? (postProcessing(@analyses, $postProcLib)) : (print "Test $testNum does not includes any post processing analysis.\n");
	} else {
		print "The file for post processing does not exists for Test $testNum.\n";
	}
	
	#~ if($observeFlag) {
		#~ if(!observableExists($specKey)) {
			#~ $observable = createObservable($tn,$specKey);
		#~ }
	#~ }
	
	
	
	#~ my $profFile = $resultsDir."prof$tn.txt";
	#~ profile($profFile);
	#~ chdir("$testDir");
	
	#~ my $dataFile = $srcDir."data$tn.txt";
	#~ my $rawFile = $testDir."raw$tn.txt";
	#~ my $energyOut = $resultsDir."e$tn.txt";
	#~ my $tstFile = $srcDir."tst$tn.txt";
	#~ my $envStack = $srcDir."EnvironStackdata$tn.txt";
	#~ my $sysStack = $srcDir."SystemStackdata$tn.txt";
	
	#~ observables($tn,$inputFile,$dataFile,$rawFile) if($observeFlag);
	
	
	
	#moveFiles($dataFile,$tstFile,$envStack,$sysStack);
	#removeFiles() if($rmFlag);
	print "******END OF TEST ".$testNum."******\n";
}

END {
	removeFiles() if($rmFlag);
}

sub postProcessing
{
	my $lib = pop(@_);
	my (@analyses) = @_;
	my %opsHash;
	
	foreach my $a(@analyses) {
		my @blockOps;
		open LIB, "<$lib" || die "Could not open file.\n";
		while(<LIB>) {
			if(/^\[$a\]/i) {
				while(<LIB>) {
					chomp;
					last if($_ eq "");
					next if(/^#/);
					push @blockOps, $_;
				}
			}
		}
		close LIB;
		
		my @commands = keyValueParser(\@blockOps, $a);		
		$opsHash{$a} = join(":", @commands);   #Grep cmd*:Diff cmd*
	}
	
	my @metaLang = ("Grep", "Make", "Call", "Execute", "Diff");
	my @recognizeCommands = ("Grep", "Make", "Diff");
	my %opsCount;  #contains @analyses count
	my %opsStack;

	foreach my $a (@analyses) {
		$opsCount{$a} = 0;
	}
	
	foreach my $a (keys %opsHash) {
		my @comm = split(/:/, "$opsHash{$a}");
		my @args;
		my @tmpLang = @metaLang;
		my $search;
		
		foreach my $key (@tmpLang) {
			$search = shift @tmpLang;			
			if(@args = grep {/$search/} @comm) {
				if(grep {/$search/} @recognizeCommands) {
					foreach my $arg (@args) {
						$arg =~ s/(\w+)/\l$1/;
						eval("hook$search(\$arg);");
					}
				} else {
					#command not recognized
				}
			}
		}
		
		$opsCount{$a}++;
	}
	
	print "Post processing has been completed.\n";
}

sub hookDiff
{
	my $rCode = system(@_);
	die "Error diff.\n" if($rCode);
	
	print "Diff successful.\n";
	
}

sub hookMake
{
	chdir($srcDir);
	my $rCode = system(@_);
	die "Error make.\n" if($rCode);
	rename "observe", "observe-$specKey";
	chdir($testDir);
	
	print "Make successful.\n";
	
}

sub hookGrep
{
	my $rCode = system(@_);
	die "Error grep.\n" if($rCode);
	
	print "Grep successful.\n";
}

sub keyValueParser
{
	my ($opsRef, $analysis) = @_;
	my @tmpKV;
	my $keyword = "Let";
	my %tmpHash;
	my @commands;
	my %varHash;
	my @varArray = ("\$srcDir", "\$inputsDir", "\$resultsDir", "\$oraclesDir", "\$testNum");
	
	foreach (@varArray) {
		$varHash{$_} = eval($_);
	}
	
	if(@tmpKV = grep {/^$keyword/} @{$opsRef}) {
		@{$opsRef} = grep {!/^$keyword/} @{$opsRef};
		grep {s/(^$keyword\s+)//g} @tmpKV;
		foreach my $keyval (@tmpKV) {
			my @t = split(/ = /, $keyval);
			$tmpHash{$t[0]} = $t[1];
		}
		
		foreach my $comm (@{$opsRef}) {
			grep {s/(\$\w+)/$tmpHash{$1}/g} $comm;
			grep {s/(\$\w+)/$varHash{$1}/g} $comm;
			grep {s/\/\s+/\//g} $comm;
			push @commands, $comm;
		}
	} else {
		die "Missing values for \"$analysis\" substitutions.\nVerify the post processing library file.\n";
	}
		
	return @commands;
}

sub extractAnalyses
{
	my ($postProc) = @_;
	my @analyses;

	open FILE, "<$postProc" || die "Could not open file.\n";
	while(<FILE>) {
		chomp;
		last if($_ eq "");
		push(@analyses, $_);
	}
	close FILE;
	
	return @analyses;
}


sub removeFiles
{
	my @files = ("Makefile*", "main*", "observe", "freeSystem*", "observe.*", "input.*", "raw$testNum.txt", "gmon.out");

	chdir($srcDir);
	system("rm @files >& /dev/null");
	chdir($testDir);
	system("rm @files >& /dev/null");
	
	print "All temporary files were removed.\n";
}

#~ sub moveFiles
#~ {
	#~ foreach (@_) { system("mv $_ ".$oraclesDir) if(validateFile("$_")); }
	
	#~ chdir("$srcDir");
	#~ foreach (@_) { system("mv $_ ".$oraclesDir) if(validateFile("$_")); }
	#~ chdir("$testDir");
	
	#~ print "All files were succesfully stored.\n";
#~ }

sub runAllTests
{
	my ($start) = @_;
	print "Preparing to run all tests starting from Test $start...\n";

	my @testsList = split(/ /,getAvailableTests());
	my @nonFunctionalTests = (4,24,60,104,105,106,107,108,124,125,141,160);
	
	for (my $i=0;$i<=$#testsList;$i++) {
		next if ($testsList[$i] eq "");
		next if ($testsList[$i]<$start);
		next if(grep {$_ eq $testsList[$i]}@nonFunctionalTests);
		testSuite($testsList[$i]);
		undef $executable;
		undef $testNum;
	}
	print "******ALL TESTS COMPLETED SUCCESSFULLY******\n";
}

