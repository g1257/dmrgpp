#! /usr/bin/perl

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
use warnings;
use Getopt::Long;
use Cwd 'abs_path';
use File::Basename;

my ($testNum) = (undef); 
my ($all,$rmFlag,$update,$verbose,$noModel,$force) = (0,0,0,0,0,0);
my ($testDir, $srcDir);
my $PATH = $testDir = $srcDir = abs_path($0);
my $filename = basename($PATH);	
$testDir =~ s/$filename$//;
$srcDir =~ s/TestSuite.*/src\//;
my $oraclesDir = $testDir."oracles/";	
my $resultsDir = $testDir."results/";
my $inputsDir = $testDir."inputs/";
my $executable = "";
my $hashTable = $testDir."hashTable.txt";
my %dmrgHash;
my %observeHash;

eval {
	my $err = GetOptions("n=i" => \$testNum, "a" => \$all, "u" => \$update, 
	"r" => \$rmFlag, "v" => \$verbose, "m" => \$noModel, "f" => \$force);
	die "Unknown command line options$!" if(!$err);
		
	&startProgram();
};
if($@) {
	die "TestSuite aborted -> $@";
}

sub startProgram{
	print "*******VALIDATION PROCESS*******\n";
	validateDirectory($srcDir);
	validateDirectory($testDir);
	validateFile($hashTable);
	validateDirectory($resultsDir);
	print "Validation process is complete.\n";
	print "*******END OF VALIDATION*******\n";
	
	if($update) {
		updateHashTables();
	} else {		
		selectTest() if(!($all) && !defined($testNum));

		loadHashTables() if(!$noModel);
		
		if($all) {
			runAllTests(0);
		} elsif($testNum < 0) {
			runAllTests(-$testNum);
		} else {
			testSuite();
		}
		
		saveHashTables() if(!$noModel);
	}
}

sub updateHashTables
{
	print "The update action for the hash table is incomplete.\n";
}

sub validateDirectory	
{
	my ($dir) = @_;
	
	if(-d $dir) {		
		return 1;
	} elsif($dir eq $resultsDir){
		mkdir($dir) || die "Making directory $dir: $!";
		$dir = basename($dir);
		print "Directory created: $dir/\n";
		return 1;
	}
	
	die "$dir: $!";
}

sub validateTest
{
	my ($available) = @_;
	my @found;
	
	if($testNum =~ /\d/) {
		my $searchNum = abs($testNum);
		@found = grep(/$searchNum/, split(/ /,$available));
	}
	
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
		$testNum = 0 if(!$testNum);
		last if(validateTest($available));
	}
}

sub getAvailableTests
{
	my $available = "";
	my $descriptionFile = $inputsDir."descriptions.txt";
	
	open(FILE,$descriptionFile) || die "Opening $descriptionFile: $!";	
	while(<FILE>) {
		last if(/^\#TAGEND/);
		next if(/^\#TAGSTART/);
		if(/(^\d+)\)/) {
			$available .= "$1 ";	
		}
		print if($verbose);	
	}
	close(FILE) || die "Closing $descriptionFile: $!";	

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
		open(FILE, ">$file") || die "Creating $file: $!";
		close(FILE) || die "Closing $file: $!";
		$file = basename($file);
		print "File created: $file\n";
		return 1;
	}
	
	die "$file: $!";
} 

sub loadHashTables
{
	open(FILE, "<$hashTable") || die "Opening $hashTable: $!";
	while(<FILE>) {
		chomp;
		if($_ eq "[dmrg]") {
			while(<FILE>) {
				chomp;
				last if($_ eq "");
				my ($key, $values) = split(/ : /);
				my @values = split(/, /, $values);
				$dmrgHash{$key} = [@values];
			}
		} elsif($_ eq "[observe]") {
			while(<FILE>) {
				chomp;
				last if($_ eq "");
				my ($key, $values) = split(/ : /);
				my @values = split(/, /, $values);
				$observeHash{$key} = [@values];
			}
		}
	}
	close(FILE) || die "Closing $hashTable: $!";
}

sub findKey
{
	my ($refHash, $key) = @_;
	
	foreach my $k (keys %{$refHash}) {
		return 1 if($key eq $k);
	}
	
	return 0;
}

sub findTn
{
	my ($refHash, $key) = @_;
	
	foreach my $k (keys %{$refHash}) {
		if($key eq $k) {
			if(grep $_ eq $testNum, @{${$refHash}{$key}}) {
				print "Hash table is up-to-date.\n";
				return 1;
			} 
		}
	}
	
	return 0;
}

sub addTn
{
	my ($refHash, $key) = @_;
	
	print "Updating hash table...\n";
	
	foreach my $k (keys %{$refHash}) {
		if($key eq $k) {
			push @{${$refHash}{$k}}, $testNum;
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
	my ($refHash, $key) = @_;
	
	print "Updating hash table...\n";
	${$refHash}{$key} = [$testNum];
}

#~ sub Keys
#~ {
	#~ my (@delKeys) = @_;
	
	#~ print "Keys removed:\n";
	#~ foreach my $key (@delKeys) {
		#~ delete $hash{$key};
		#~ print "$key\n";
	#~ }
#~ }

sub saveHashTables
{
	open(FILE, ">$hashTable") || die "Opening $hashTable: $!";
	print FILE "[dmrg]\n";
	foreach my $key (sort keys %dmrgHash) {
		print FILE "$key : ".join(', ', sort @{$dmrgHash{$key}})."\n";
	}
	print FILE "\n";
	print FILE "[observe]\n";
	foreach my $key (sort keys %observeHash) {
		print FILE "$key : ".join(', ', sort @{$observeHash{$key}})."\n";
	}
	close(FILE) || "Closing $hashTable: $!";
}

sub createExecutable
{
	my ($specFile,$specKey,$configFile) = @_;
	my $arg1 = "./$configFile < $specFile >& /dev/null";
	my $arg2 = "make dmrg -f Makefile >& /dev/null";
	
	grep {s/<.*//} $arg1 if($noModel);
	grep {s/>.*//} $arg1 if($verbose);
	grep {s/>.*//} $arg2 if($verbose);
	
	my $err = chdir($srcDir);
	die "Changing directory to $srcDir: $!" if(!$err);
	my $err = system($arg1);
	die "Configuration error using $configFile with $specFile: $!" if($err);
	print "Configuration of Test $testNum was successful.\n";
	print "Creating executable for Test $testNum...\n";
	my $err = system($arg2);
	die "Make command for dmrg: $!" if($err);
	if(!$noModel) {
		my $oldExe = $executable;
		$executable = $executable."-".$specKey;
		my $err = rename($oldExe, $executable);
		die "Renaming $oldExe to $executable: $!" if(!$err);
	}
	my $err = chdir($testDir);
	die "Changing directory to $testDir: $!" if(!$err);
	
	print "Executable was succesfully created.\n";
	
	return $srcDir.$executable;	
}

sub getSpecFileAndKey
{
	my $tempNum = $testNum;
	$tempNum -= 100 if($testNum >= 100);
	my $specFile = $inputsDir."model$tempNum.spec";
	
	my $specKey = substr(`md5sum $specFile`,0,10);
	
	return $specFile, $specKey;
}

sub testSuite
{
	my $tempNum = $testNum;
	$tempNum -= 100 if($testNum >= 100);
	my $procFile = $inputsDir."processing$tempNum.txt";
	my $procLib = $inputsDir."processingLibrary.txt";
	
	print "*******START OF TEST $testNum*******\n";
	
	if(validateFile($procLib)) {
		if(validateFile($procFile)) {
			my @analyses = extractAnalyses($procFile) ;
			(@analyses) ? (processing(@analyses, $procLib)) : (print "Test $testNum does not includes any processing analyses.\n");
		} 
	}
	
	removeFiles() if($rmFlag && !$update);
	print "*******END OF TEST ".$testNum."*******\n";
}

sub processing
{
	my $lib = pop(@_);
	my (@analyses) = @_;
	my %procHash;
	my @commands;
	
	foreach my $analysis(@analyses) {
		my @operations;
		open(LIB, "<$lib") || die "Opening $lib: $!";
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
		close(LIB) || die "Closing $lib: $!";
		
		if(!@operations) {
			print "No commands are described in Test $testNum for [$analysis] analysis.\n";
		} else {
			my @commands = keyValueParser(\@operations);		
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
			
			die "No runable command in library for [$analysis] analysis.\n" if(!@commands);
			commandsInterpreter(@commands, $analysis);
			$procCount{$analysis}++;
			delete $tmpHash{$analysis};
		}

		die "Unresolved processing dependencies in $lib.\n" if($prevCount == keys %tmpHash);
	}
	
	print "Processing completion successful.\n";
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
			if($@) {
				my $subr = caller(0)[3];
				die "Subroutine $subr: $@";
			}
		} else {
			die "Unknown command in library analysis [$analysis]: $func\n";
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
	if($@) {
		my $subr = caller(0)[3];
		die "Subroutine $subr: $@";
	}
}

sub runDmrg
{
	my ($inputFile) = @_;
	my ($specFile, $specKey) = getSpecFileAndKey();
	my $configFile = "configure.pl";
	my $dataFile = "data$testNum.txt";
	$executable = "dmrg";
	
	if(!findKey(\%dmrgHash, $specKey) || $force || $noModel) {
		$executable = createExecutable($specFile,$specKey,$configFile);
		addKey(\%dmrgHash, $specKey);
	} else {
		print "Retrieving existing executable...\n";
		$executable = $srcDir.$executable."-".$specKey;
		validateFile($executable);
	
		if(!findTn(\%dmrgHash, $specKey) ){
			addTn(\%dmrgHash, $specKey);
		}
	}
	
	my $arg = "$executable $inputFile >& /dev/null";
	grep {s/&//} $arg if($verbose);
	
	print "Please wait while the test is ran...\n";
	my $err = chdir($srcDir);
	die "Changing directory to $srcDir: $!" if(!$err);
	my $err = system($arg);
	die "Running test using $executable with $inputFile: $!" if($err);
	my $err = chdir($testDir);
	die "Changing directory to $testDir: $!" if(!$err);
	print "The run has been completed.\n";
	
	#move("$srcDir$dataFile", "$resultsDir$dataFile");
}

sub runObserve
{
	my ($input, $raw) = @_;
	my ($specFile, $specKey) = getSpecFileAndKey();
	my $configFile = "configure.pl";
	my $observable = "observe";
	
	my $err = chdir($srcDir);
	die "Changing directory to $srcDir: $!" if(!$err);
	if(!validateFile("Makefile")) {
		my $arg1 = "./$configFile < $specFile >& /dev/null";
		grep {s/<.*//} $arg1 if($noModel);
		grep {s/>.*//} $arg1 if($verbose);
		my $err = system($arg1);
		die "Configuration error using $configFile with $specFile: $!" if($err);
		print "Configuration of Test $testNum was successful...\n";
	}
	
	my $arg2 = "make observe -f Makefile >& /dev/null";
	grep {s/>.*//} $arg2 if($verbose);
	
	print "Creating executable for Test $testNum...\n";
	my $err = system($arg2);
	die "Make command for observables: $!" if($err);
	if(!$noModel) {
		my $oldObs = $observable;
		$observable = $observable."-".$specKey;
		my $err = rename($oldObs, $observable);
		die "Renaming $oldObs to $observable: $!" if(!$err);
	}
	
	print "Observable was succesfully created.\n";
	my $arg = "./$observable $input > $raw";
	print "Please wait while the observable is run...\n";
	my $err = system($arg);
	die "Running test using $observable with $input: $!" if($err);
	my $err = chdir($testDir);
	die "Changing directory to $testDir: $!" if(!$err);
	print "The run has been completed.\n";
}

sub hookGprof
{
	my ($analysis, $arg) = @_;
	
	my $err = chdir($srcDir);
	die "Changing directory to $srcDir: $!" if(!$err);
	eval("system(\"gprof $arg\");");
	if($@) {
		my $subr = caller(0)[3];
		die "Subroutine $subr: $@";
	}
	my $err = chdir($testDir);
	die "Changing directory to $testDir: $!" if(!$err);
	
	print "$analysis:Gprof was successful.\n";
	
}

sub hookDiff
{
	my ($analysis, $arg) = @_;
	
	eval("system(\"diff $arg\");");
	if($@) {
		my $subr = caller(0)[3];
		die "Subroutine $subr: $@";
	}

	print "$analysis:Diff was successful.\n";
	
}

sub hookMake
{
	my ($analysis, $arg) = @_;
	
	my $err = chdir($srcDir);
	die "Changing directory to $srcDir: $!" if(!$err);
	eval("system(\"make $arg\");");
	if($@) {
		my $subr = caller(0)[3];
		die "Subroutine $subr: $@";
	}
	my $err = chdir($testDir);
	die "Changing directory to $testDir: $!" if(!$err);
	
	print "$analysis:Make was successful.\n";
}

sub hookGrep
{
	my ($analysis, $arg) = @_;
	
	eval("system(\"grep $arg\");");
	if($@) {
		my $subr = caller(0)[3];
		die "Subroutine $subr: $@";
	}
	
	print "$analysis:Grep was successful.\n";
}

sub extractOperatorC
{
	my ($raw,$cOut) = @_;
	my $line;
	my $opC;
	
	open(INFILE,"<$raw") || die "Opening $raw: $!";
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
	close(INFILE) || die "Closing $raw: $!";
	
	open (OUTFILE, ">$cOut") || die "Opening $cOut: $!";
	print OUTFILE $opC;
	close (OUTFILE) || die "Closing $cOut: $!";
	print "OperatorC extraction was succesful!\n";
}

sub extractOperatorN
{
	my ($raw,$nOut) = @_;
	my $line;
	my $opN;
	open(INFILE,"<$raw") || die "Opening $raw: $!";
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
	close(INFILE) || die "Closing $raw: $!";
	
	open (OUTFILE, ">$nOut") || die "Opening $nOut: $!";
	print OUTFILE $opN;
	close (OUTFILE) || die "Closing $nOut: $!";
	print "OperatorN extraction was succesful!\n";
}


sub extractOperatorSz
{
	my ($raw,$sOut) = @_;
	my $line;
	my $opS;
	open(INFILE,"<$raw") || die "Opening $raw: $!";
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
	close(INFILE) || die "Closing $raw: $!";
	
	open (OUTFILE, ">$sOut") || die "Opening $sOut: $!";
	print OUTFILE $opS;
	close (OUTFILE) || die "Closing $sOut: $!";
	print "OperatorSz extraction was succesful!\n";
}

sub keyValueParser
{
	my ($opsRef) = @_;
	my @tmpKV;
	my $keyword = "Let";
	my %tmpHash;
	my @commands;
	my %varHash;
	my @varArray = ("\$executable", "\$srcDir", "\$inputsDir", "\$resultsDir", "\$oraclesDir", "\$testNum");
	my @nonSubsArray = ("\$executable");
	
	foreach (@varArray) {
		my $tmp = "\\$_";
		if(grep {/$tmp/} @nonSubsArray) {
			$varHash{$_} = $_;
		} else {
			$varHash{$_} = eval($_);
		}
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
			grep {s/(\$\w+)/$tmpHash{$1}/g} $comm;
			grep {s/(\$\w+)(\s+)([^<>])/$varHash{$1}$3/g} $comm;
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

	open(FILE, "<$procFile") || die "Opening $procFile: $!";
	while(<FILE>) {
		chomp;
		last if($_ eq "");
		next if(/^#/);
		push(@analyses, $_);
	}
	close (FILE) || die "Closing $procFile: $!";
	
	return @analyses;
}


sub removeFiles
{
	my @files = ("Makefile*", "main*", "observe", "freeSystem*", "observe.*", "input.*", "raw$testNum.txt", "gmon.out");

	my $err = chdir($srcDir);
	die "Changing directory to $srcDir: $!" if(!$err);
	my $err = system("rm @files >& /dev/null");
	die "Removing files: $!" if($err);
	my $err = chdir($testDir);
	die "Changing directory to $testDir: $!" if(!$err);
	my $err = system("rm @files >& /dev/null");
	die "Removing files: $!" if($err);
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
		$testNum = $testsList[$i];
		testSuite();
		$executable = "";
		removeFiles() if($rmFlag);
	}
	print "******ALL TESTS COMPLETED SUCCESSFULLY******\n";
}

