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
#use File::Copy;

my ($testNum,$all) = (undef,undef); 
my ($rmFlag,$update,$verbose,$noModel,$force) = (0,0,0,0,0);

GetOptions("n=i" => \$testNum, "a" => \$all, "u" => \$update, 
"r" => \$rmFlag, "v" => \$verbose, "m" => \$noModel, "f" => \$force);

my ($testDir, $srcDir);
my $PATH = $testDir = $srcDir = abs_path($0);
my $filename = basename($PATH);	

$testDir =~ s/$filename$//;
$srcDir =~ s/TestSuite.*/src\//;
my $oraclesDir = $testDir."oracles/";	
my $resultsDir = $testDir."results/";
my $inputsDir = $testDir."inputs/";

my $executable ;
my $hashTable = $testDir."hashTable.txt";
my %hash;

&startProgram();

sub startProgram{
	if($update) {
		updateHashTable();
	} elsif(validateDirectory($srcDir) && validateDirectory($testDir) && validateDirectory($resultsDir) && validateFile($hashTable)) {
		
		selectTest() if(!($all) && !defined($testNum));

		loadHashTable() if(!$noModel);

		if($testNum < 0) {
			runAllTests(-$testNum);
		} elsif($all) {
			runAllTests(0);
		} else {
			testSuite();
		}
		
		saveHashTable() if(!$noModel);
			
	} else {
		die "Missing source and/or test suite directories.\n";
	}
}

sub updateHashTable
{
	print "The update action for the hash table is incomplete.\n";
}

sub validateDirectory	
{
	my ($dir) = @_;
	
	if(-d $dir) {		
		return 1;
	} elsif($dir eq $resultsDir){
		my $rCode = mkdir $dir;
		if($rCode) {
			print "Directory created: $dir\n";
			return 1;
		}
	}
	
	return 0;
}

sub validateTest
{
	my ($available) = @_;
	
	my $searchNum = abs($testNum);
	my @found = grep(/$searchNum/, split(/ /,$available));	
	
	if(!@found) {	
		print "\nError: An incorrect test was selected! Try again.\n\n";
		return 0;
	}

	return 1;
}

sub selectTest
{
	my $available = getAvailableTests();
	
	while() {
		print "Type the number of the test you want to run.\n";
		print "Available tests: $available\n";
		print "Default is 0 (press ENTER): ";
		chomp($testNum = <STDIN>);
		$testNum = 0 if($testNum eq "");
		last if(validateTest($available));
	}
}

sub getAvailableTests
{
	my $available = "";
	my $descriptionFile = $inputsDir."descriptions.txt";
	
	open(FILE,$descriptionFile) || die "Error opening $descriptionFile: $!\n";	
	while(<FILE>) {
		last if(/^\#TAGEND/);
		next if(/^\#TAGSTART/);
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
		my $rCode = open(FILE, ">$hashTable") || die "Error creating $hashTable: $!\n";
		close(FILE);
		if($rCode) {
			print "File created: $hashTable\n";
			return 1;
		}
	}
	
	return 0;
} 

sub loadHashTable
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

sub saveHashTable
{
	open FILE, ">$hashTable" || die "Could not open file: $hashTable\n";
	foreach my $key (sort keys %hash) {
		print FILE "$key : ".join(', ', sort @{$hash{$key}})."\n";
	}
	close FILE;
}

sub createExecutable
{
	my ($specFile,$specKey,$configFile) = @_;
	my $arg1 = "./$configFile < $specFile >& /dev/null";
	my $arg2 = "make -f Makefile >& /dev/null";
	
	grep {s/<.*//} $arg1 if($noModel);
	grep {s/>.*//} $arg1 if($verbose);
	grep {s/>.*//} $arg2 if($verbose);
	
	chdir($srcDir);
	system($arg1);
	print "Configuration of Test $testNum was successful...\n";
	print "Creating executable for Test $testNum...\n";
	system($arg2);
	if(!$noModel) {
		my $oldExe = $executable;
		$executable = $executable."-".$specKey;
		rename($oldExe, $executable);
	}
	chdir($testDir);
	
	print "Executable was succesfully created.\n";
	
	return $srcDir.$executable;	
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

#~ sub testSuite
#~ {
	#~ my $temp = $testNum;
	#~ $temp -= 100 if($testNum >= 100);
	#~ my $specFile = $inputsDir."model$temp.spec";
	#~ my $configFile = "configure.pl";
	#~ my $inputFile = $inputsDir."input$testNum.inp";
	#~ my $procFile = $inputsDir."processing$testNum.txt";
	#~ my $procLib = $inputsDir."processingLibrary.txt";
	
	#~ $specKey = substr(`md5sum $specFile`,0,10);
	
	#~ if(defined($executable)){
		#~ print "Using executable supplied.\n";
		#~ $executable = $srcDir.$executable;
		#~ die "Missing executable file: $executable" if(!validateFile($executable));
	#~ }elsif($noModel) {
		#~ $executable = createManualExecutable($configFile);
	#~ }elsif(!findKey($specKey) || $force) {
		#~ $executable = createAutoExecutable($specFile,$specKey,$configFile);
		#~ addKey($specKey);
	#~ } else {
		#~ print "Retrieving existing executable...\n";
		#~ $executable = $srcDir."dmrg-$specKey";
		#~ validateFile($executable) || die "Error retrieving $executable: $!\n";
	
		#~ if(!findTn($specKey) ){
			#~ addTn($specKey);
		#~ }
	#~ }
	
	#~ #runTest($inputFile);
	#~ if(validateFile($procLib)) {
		#~ if(validateFile($procFile)) {
			#~ my @analyses = extractAnalyses($procFile) ;
			#~ (@analyses) ? (processing(@analyses, $procLib)) : (print "Test $testNum does not includes any processing analyses.\n");
		#~ } else {
			#~ print "The file for processing does not exists for Test $testNum.\n";
		#~ }
	#~ } else {
		#~ print "The library for processing does not exists.\n";
	#~ }
	
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
	
	
	
	#~ #moveFiles($dataFile,$tstFile,$envStack,$sysStack);
	#~ #removeFiles() if($rmFlag);
	#~ print "******END OF TEST ".$testNum."******\n";
#~ }

sub getSpecFileAndKey
{
	my $temp = $testNum;
	$temp -= 100 if($testNum >= 100);
	my $specFile = $inputsDir."model$temp.spec";
	
	my $specKey = substr(`md5sum $specFile`,0,10);
	
	return $specFile, $specKey;
}

sub testSuite
{
	my $procFile = $inputsDir."processing$testNum.txt";
	my $procLib = $inputsDir."processingLibrary.txt";
	
	if(validateFile($procLib)) {
		if(validateFile($procFile)) {
			my @analyses = extractAnalyses($procFile) ;
			(@analyses) ? (processing(@analyses, $procLib)) : (print "Test $testNum does not includes any processing analyses.\n");
		} else {
			print "The file for processing Test $testNum does not exists.\n";
		}
	} else {
		print "The library for processing does not exists.\n";
	}

	print "******END OF TEST ".$testNum."******\n";
}

END {
	removeFiles() if($rmFlag && !$update);
}

sub processing
{
	my $lib = pop(@_);
	my (@analyses) = @_;
	my %procHash;
	my @commands;
	
	foreach my $analysis(@analyses) {
		my @operations;
		open LIB, "<$lib" || die "Could not open file.\n";
		while(<LIB>) {
			if(/^\[$analysis\]/i) {
				while(<LIB>) {
					chomp;
					last if($_ eq "");
					next if(/^#/);
					push @operations, $_;
				}
			}
		}
		close LIB;
		
		if(!@operations) {
			print "No commands are described in Test $testNum for [$analysis] analysis.\n";
		} else {
			my @commands = keyValueParser(\@operations, $analysis);		
			$procHash{$analysis} = join(":", @commands);
		}
	}
	
	die "No analysis was found for Test $testNum.\n" if(!keys %procHash);
	
	my %procCount;
	my %tmpHash = %procHash;
	my @dependencies;
	my $depFlag;
	my @depKeys;
	my $keyword = "CallOnce";
	my $prevCount = 0;

	foreach my $analysis (keys %tmpHash) {
		$procCount{$analysis} = 0;
	}
	
	while($prevCount = keys %tmpHash) {
		foreach my $analysis( keys %tmpHash) {
			$depFlag = 0;
			@commands = split(/:/, $tmpHash{$analysis});
			if(@dependencies = grep{/$keyword/} @commands) {
				@depKeys = grep{s/$keyword|\s+//g} @dependencies;
				foreach my $a(@depKeys) {
					if($procCount{$a} eq 0) {
						$depFlag = 1;
					}
				}

				if($depFlag) {
					next;
				} else {
					@commands = grep {!/$keyword/} @commands;
				}
			}
			
			die "Error due to no runable command in library analysis [$analysis].\n" if(!@commands);
			commandsInterpreter(@commands, $analysis);
			$procCount{$analysis}++;
			delete $tmpHash{$analysis};
		}

		die "Error with processing dependencies.\n" if($prevCount == keys %tmpHash);
	}
	
	print "Processing has been completed.\n";
}

sub commandsInterpreter
{
	my $analysis = pop(@_);
	my (@commands) = @_;
	my @metaLang = ("Grep", "Make", "Execute", "Gprof", "Diff");
	my @arrangeCommands;
	
	foreach my $word(@metaLang) {
		if(my @tmpComm = grep {/$word/} @commands) {
			foreach my $comm (@tmpComm) {
				push @arrangeCommands, $comm;
			}
		}
	}
	
	foreach my $arg (@arrangeCommands) {
		my @tmpFunc = $arg =~ /^(\w+)/;
		my $func = $tmpFunc[0];
		$arg =~ s/^\s*(\w+)\s*//;
		
		if(grep {/$func/} @metaLang) {
			eval("hook$func(\$analysis,\$arg);");
		} else {
			die "Error with unknown command in library analysis [$analysis]: $func\n";
		}
	}
}

sub hookExecute
{
	my ($analysis, $arg) = @_;
	
	$arg =~ s/(\()/$1"/g;
	$arg =~ s/(\s*)(,)(\s*)/"$2"/g;
	$arg =~ s/(\))/"$1/g;
	
	eval("$arg;");
	warn $@ if $@;
}

sub runDmrg
{
	my ($inputFile) = @_;
	my ($specFile, $specKey) = getSpecFileAndKey();
	my $configFile = "configure.pl";
	my $dataFile = "data$testNum.txt";
	$executable = "dmrg";
	
	if(!findKey($specKey) || $force || $noModel) {
		$executable = createExecutable($specFile,$specKey,$configFile);
		addKey($specKey);
	} else {
		print "Retrieving existing executable...\n";
		$executable = $srcDir.$executable."-".$specKey;
		validateFile($executable) || die "Error retrieving $executable: $!\n";
	
		if(!findTn($specKey) ){
			addTn($specKey);
		}
	}
	
	runTest($inputFile);
	system("mv $srcDir$dataFile $resultsDir$dataFile");
}

sub runObserve
{
	my ($input, $raw) = @_;
	
	system("../src/observe $input > $raw");
	print "run observe\n";
	
}

sub hookGprof
{
	my ($analysis, $arg) = @_;
	
	my $rCode = system("gprof $arg");
	die "Error gprof.\n" if($rCode);
	
	print "Gprof was successful.\n";
	
}

sub hookDiff
{
	my ($analysis, $arg) = @_;
	
	my $rCode = system("diff $arg");
	die "Error diff.\n" if($rCode);
	
	print "Diff was successful.\n";
	
}

sub hookMake
{
	my ($analysis, $arg) = @_;
	
	chdir($srcDir);
	my $rCode = system("make $arg");
	die "Error make.\n" if($rCode);
	chdir($testDir);
	
	print "Make was successful.\n";
}

sub hookGrep
{
	my ($analysis, $arg) = @_;
	
	my $rCode = system("grep $arg");
	die "Error grep.\n" if($rCode);
	
	print "Grep was successful.\n";
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
	}	
	
	foreach my $comm (@{$opsRef}) {
		if(@tmpKV) {
			#me falta validar los $key/values
			grep {s/(\$\w+)/$tmpHash{$1}/g} $comm;
			grep {s/(\$\w+)\s+/$varHash{$1}/g} $comm;
			grep {s/(\$\w+)/$varHash{$1}/g} $comm;
		}
		push @commands, $comm;
	}
		
	return @commands;
}

sub extractAnalyses
{
	my ($procFile) = @_;
	my @analyses;

	open FILE, "<$procFile" || die "Could not open file.\n";
	while(<FILE>) {
		chomp;
		last if($_ eq "");
		next if(/^#/);
		push(@analyses, $_);
	}
	close FILE;
	
	die "No analyses were found in $procFile for Test $testNum\n" if(!@analyses);
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

