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

#Declarations of Perl modules required
use strict;
use warnings;
use Getopt::Long;
use Cwd 'abs_path';

#Global variables and command line options
my ($testNum,$lastTest) = ("",-1); 
my ($all,$rmFlag,$update,$verbose,$noModel,$force) = (0,0,0,0,0,0);
my ($testDir, $srcDir);
my $PATH = $testDir = $srcDir = abs_path($0);
chomp(my $filename = `basename $0`);
$testDir =~ s/$filename$//;
$srcDir =~ s/TestSuite.*/src\//;
my $oraclesDir = $testDir."oracles/";	
my $resultsDir = $testDir."results/";
my $inputsDir = $testDir."inputs/";
my $executable = "";
my $hashTable = $testDir."hashTable.txt";
my %dmrgHash = ();
my %observeHash = ();
my $globalRunFlag = 0;

#Signals for Control+C escape and warnings
$SIG{INT} = \&exit_handler;
$SIG{__WARN__} = sub {die @_};

#Handler for Control+C escape signal
sub exit_handler
{
	print "\nTestSuite aborted -> Manual cancellation\n";
	removeFiles();
	kill 9, $$;
}

eval {
	#Get command line options
	die $! if(!GetOptions("n=i" => \$testNum, "l=i" => \$lastTest, "all" => \$all, "u" => \$update, 
	"r" => \$rmFlag, "v" => \$verbose, "m" => \$noModel, "f" => \$force));
	
	#Activate testsuite program
	$globalRunFlag = 1;
	
	while($globalRunFlag) {
		#Starts main program
		startProgram();
	}
};
if($@) {
	print "\nTestSuite aborted -> $@";
	removeFiles();
	kill 9, $$;
}

sub startProgram{
	print "*******INITIAL PROCESSING*******\n";
	die "$!" if(!validateDirectory($srcDir));
	die "$!" if(!validateDirectory($testDir));
	die "$!" if(!validateFile($hashTable));
	die "$!" if(!validateDirectory($resultsDir));
	
	if($update) {
		updateHashTables();
	} else {		
		selectTest() if(!($all) && ($testNum eq ""));

		if($all) {
			runAllTests(0);
		} elsif($testNum < 0) {
			runAllTests(-$testNum);
		} else {
			testSuite();
		}
		
		print "*******FINAL PROCESSING*******\n";
	}
	
	$globalRunFlag = 0;
}

sub updateHashTables
{
	my $err = chdir($srcDir);
	die "Changing directory to $srcDir: $!" if(!$err);
	my @dmrgFiles = glob("dmrg-*");
	my @observeFiles = glob("observe-*");
	my @dmrgKeys = grep{s/^(dmrg-)//g} @dmrgFiles;
	my @observeKeys = grep{s/^(observe-)//g} @observeFiles;	
	$err = chdir($testDir);
	die "Changing directory to $testDir: $!" if(!$err);
	
	print "*******HASH TABLE UPDATE*******\n";
	loadHashTable();
	my $deleteFlag;
	
	foreach my $key (keys %dmrgHash) {
		$deleteFlag = 1;
		
		if(grep{$_ eq $key} @dmrgKeys) {
			$deleteFlag = 0;
		}
		
		delete $dmrgHash{$key} if($deleteFlag);
	}
	
	foreach my $key(@dmrgKeys) {
		if(!findKey(\%dmrgHash, $key)) {
			addKey(\%dmrgHash, $key);
		}
	}
	
	foreach my $key (keys %observeHash) {
		$deleteFlag = 1;
		
		if(grep{$_ eq $key} @observeKeys) {
			$deleteFlag = 0;
		}
		
		delete $observeHash{$key} if($deleteFlag);
	}
	
	foreach my $key(@observeKeys) {
		if(!findKey(\%observeHash, $key)) {
			addKey(\%observeHash, $key);
		}
	}
	
	saveHashTable();
	print "*******END OF UPDATE*******\n";
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

sub validateTest
{
	my ($available) = @_;
	
	if($testNum =~ /(^-?\d+$)/) {
		my $searchNum = abs($testNum);
		my @found = grep(/$searchNum/, split(/ /,$available));
		if(@found) {
			if($testNum eq "-0") {
				$testNum = 0;
				$all = 1;
			}
			return 1;
		}
	}	
	print "\n<Error>: An incorrect test was selected! Try again.\n\n";
	return 0;
}

sub validateDirectory	
{
	my ($dir) = @_;
	
	if(-d $dir) {		
		return 1;
	} elsif($dir eq $resultsDir){
		mkdir($dir) || die "Making directory $dir: $!";
		print "Directory created: $dir/\n";
		return 1;
	}
	
	return 0;
}

sub validateFile	
{
	my ($file) = @_;
	
	if(-e $file) {		
		return 1;
	} elsif($file eq $hashTable) {
		open(FILE, ">$file") || die "Creating $file: $!";
		close(FILE) || die "Closing $file: $!";
		print "File created: $file\n";
		return 1;
	}
	
	return 0;
} 

sub loadHashTable
{
	open(FILE, "<$hashTable") || die "Opening $hashTable: $!";
	while(<FILE>) {
		chomp;
		if($_ eq "[dmrg]") {
			while(<FILE>) {
				chomp;
				last if($_ eq "");
				my ($key, $values) = split(/ : /);
				$values = "" if(!$values && $values ne 0);
				my @values = split(/, /, $values);
				$dmrgHash{$key} = [@values];
			}
		} elsif($_ eq "[observe]") {
			while(<FILE>) {
				chomp;
				last if($_ eq "");
				my ($key, $values) = split(/ : /);
				$values = "" if(!$values && $values ne 0);
				my @values = split(/, /, $values);
				$observeHash{$key} = [@values];
			}
		}
	}
	close(FILE) || die "Closing $hashTable: $!";
	
	print "Loading of hash table was successful.\n";
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

sub addKey
{
	my ($refHash, $key) = @_;
	
	${$refHash}{$key} = [$testNum];
}

sub addTn
{
	my ($refHash, $key) = @_;
	
	foreach my $k (keys %{$refHash}) {
		push (@{${$refHash}{$k}}, $testNum) if($key eq $k);
	}
}

sub saveHashTable
{
	open(FILE, ">$hashTable") || die "Opening $hashTable: $!";
	if(keys %dmrgHash) {
		print FILE "[dmrg]\n";
		foreach my $key (sort keys %dmrgHash) {
			print FILE "$key : ".join(', ', sort @{$dmrgHash{$key}})."\n";
		}
		print FILE "\n";
	}
	
	if(keys %observeHash) {
		print FILE "[observe]\n";
		foreach my $key (sort keys %observeHash) {
			print FILE "$key : ".join(', ', sort @{$observeHash{$key}})."\n";
		}
	}
	close(FILE) || die "Closing $hashTable: $!";
	
	print "Saving of hash table was successful.\n";
}

sub testSuite
{
	my $tempNum = $testNum;
	$tempNum -= 100 if($testNum >= 100);
	my $procFile = $inputsDir."processing$tempNum.txt";
	my $procLib = $inputsDir."processingLibrary.txt";
	
	if(validateFile($procLib)) {
		if(validateFile($procFile)) {
			print "*******START OF TEST $testNum*******\n";
			loadHashTable();
			my @analyses = extractAnalyses($procFile) ;
			(@analyses) ? (processing(@analyses, $procLib)) : (print "Test $testNum does not includes any processing analyses.\n");
			saveHashTable();
			print "*******END OF TEST ".$testNum."*******\n";
		} else {
			die "$!";
		}
	} else {
		die "$!";
	}
	
	moveFiles();
	removeFiles() if($rmFlag);
	print "All temporary files were successfully removed.\n";
}

sub runAllTests
{
	my ($start) = @_;
	my @nonFunctionalTests = (24,60,104,105,106,107,108,124,125,141,160);
	print "Preparing to run all tests starting from Test $start...\n";

	my @testsList = split(/ /,getAvailableTests());
	
	for (my $i=0;$i<=$#testsList;$i++) {
		next if ($testsList[$i] eq "");
		next if ($testsList[$i]<$start);
		next if(grep {$_ eq $testsList[$i]}@nonFunctionalTests);
		$testNum = $testsList[$i];
		testSuite();
		$executable = "";
		last if($testsList[$i] == $lastTest);
	}
	print "******ALL TESTS COMPLETED SUCCESSFULLY******\n";
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
					last if($_ eq "" || /(^\[)*(\])/);
					next if(/^#/);
					push @operations, $_;
				}
			}
		}
		close(LIB) || die "Closing $lib: $!";
		
		if(!@operations) {
			print "<Warning>: Missing descriptions for [$analysis] analysis in $lib.\n";
		} else {
			my @commands = keyValueParser(\@operations);
			if(!@commands) {
				print "<Warning>: No runnable commands for [$analysis] analysis in $lib.>\n";
			} else {
				$procHash{$analysis} = join(":", @commands);
			}
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
					if(!defined($procCount{$a})) {
						die "Unresolved dependency for analysis [$analysis]: inactive analysis [$a]. Verify the processing file for Test $testNum or processing library.\n";
					} elsif($procCount{$a} eq 0) {
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

sub commandsInterpreter
{
	my $analysis = pop(@_);
	my (@commands) = @_;
	my @metaLang = ("Grep", "Execute", "Gprof", "Diff");
	my @arrangeCommands;
	
	foreach my $word(@metaLang) {
		if(my @tmpComm = grep {/^$word/} @commands) {
			@commands = grep {!/^$word/} @commands;
			foreach my $comm (@tmpComm) {
				push @arrangeCommands, $comm;
			}
		}
	}
	
	die "Unknown commands in analysis [$analysis]: @commands\n" if(@commands);
	
	foreach my $arg (@arrangeCommands) {
		my @tmpFunc = $arg =~ /^(\w+)/;
		$arg =~ s/^\s*(\w+)\s*//;
		eval("hook$tmpFunc[0](\$analysis,\$arg);");
		if($@) {
			my $subr = (caller(0))[3];
			die "$subr: $@";
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
		my $subr = (caller(0))[3];
		die "$subr: $@";
	}
}

sub runDmrg
{
	my ($inputFile,$raw) = @_;
	my ($specFile, $specKey) = ("", "");
	($specFile, $specKey) = getSpecFileAndKey() if(!$noModel);
	my $configFile = "configure.pl";
	$executable = "dmrg";

	if(!findKey(\%dmrgHash, $specKey) || $force || $noModel) {
		$executable = createExecutables($specFile,\$specKey,$configFile, $executable);
		print "Updating hash table...\n";
		addKey(\%dmrgHash, $specKey);
	} else {
		print "Retrieving existing executable...\n";
		$executable = $srcDir.$executable."-".$specKey;
		die "$!" if(!validateFile($executable));
	
		if(!findTn(\%dmrgHash, $specKey) ){
			print "Updating hash table...\n";
			addTn(\%dmrgHash, $specKey);
		}
	}

	my $arg = "$executable $inputFile >& $raw";
	grep {s/&//} $arg if($verbose);
	
	print "Running dmrg test...\n";
	my $err = chdir($srcDir);
	die "Changing directory to $srcDir: $!" if(!$err);
	$err = system($arg);
	die "Running test using $executable with $inputFile: $!" if($err);
	$err = chdir($testDir);
	die "Changing directory to $testDir: $!" if(!$err);
	print "The dmrg run has been completed.\n";
}

sub runObserve
{
	my ($inputFile, $raw) = @_;
	my ($specFile, $specKey) = ("", "");
	($specFile, $specKey) = getSpecFileAndKey() if(!$noModel);
	my $configFile = "configure.pl";
	$executable = "observe";
	
	if(!findKey(\%observeHash, $specKey) || $force || $noModel) {
		$executable = createExecutables($specFile,\$specKey,$configFile, $executable);
		print "Updating hash table...\n";
		addKey(\%observeHash, $specKey);
	} else {
		print "Retrieving existing executable...\n";
		$executable = $srcDir.$executable."-".$specKey;
		die "$!" if(!validateFile($executable));
	
		if(!findTn(\%observeHash, $specKey) ){
			print "Updating hash table...\n";
			addTn(\%observeHash, $specKey);
		}
	}
	
	my $arg = "$executable $inputFile >& $raw";
	grep {s/&//} $arg if($verbose);
	
	print "Running observe test...\n";
	my $err = chdir($srcDir);
	die "Changing directory to $srcDir: $!" if(!$err);
	$err = system($arg);
	die "Running test using $executable with $inputFile: $!" if($err);
	$err = chdir($testDir);
	die "Changing directory to $testDir: $!" if(!$err);
	print "The observe run has been completed.\n";
}

sub createExecutables
{
	my ($specFile,$refKey,$configFile, $execType) = @_;
	my $arg1 = "./$configFile < $specFile >& /dev/null";
	my $arg2 = "make $execType -f Makefile >& /dev/null";

	grep {s/<.*//} $arg1 if($noModel);
	grep {s/>.*//} $arg1 if($verbose);
	grep {s/>.*//} $arg2 if($verbose);
	
	my $err = chdir($srcDir);
	die "Changing directory to $srcDir: $!" if(!$err);
	$err = system($arg1);
	die "Configuration error using $configFile with $specFile: $!" if($err);
	print "Configuration of $execType Test $testNum was successful.\n";
	print "Creating $execType executable for Test $testNum...\n";
	$err = system($arg2);
	die "Make command for $execType: $!" if($err);
	
	if($noModel) {
		while() {
			print "Enter a unique alphanumeric key for the $execType executable: ";
			chomp(${$refKey} = <STDIN>);
			if(${$refKey} eq "") {
				print "\n";
				next;
			}
			next if(grep{/\s/} ${$refKey});
			next if(grep{/\W+/} ${$refKey});
			last;
		}
	}
	
	my $oldExec = $executable;
	$executable= $executable."-".${$refKey};
	$err = rename($oldExec, $executable);
	die "Renaming $oldExec to $executable: $!" if(!$err);
	
	$err = chdir($testDir);
	die "Changing directory to $testDir: $!" if(!$err);
	
	print "\u$execType executable was succesfully created.\n";
	
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

sub extractOperator
{
	my ($opName, $raw,$out) = @_;
	my $line;
	my $op;
	
	open(INFILE,"<$raw") || die "Opening $raw: $!";
		while($line = <INFILE>) {
			if($line =~ /^Operator$opName/) {
				$op = $line;
				$line = <INFILE>;
				$op = $op.$line;
				my @temp1 = split(/ /,$line);

				for(my $i = 0; $i < $temp1[0]; $i++) {
					$line = <INFILE>;
					$op = $op.$line;
				}
			}
		}
	close(INFILE) || die "Closing $raw: $!";
	
	open (OUTFILE, ">$out") || die "Opening $out: $!";
	print OUTFILE $op;
	close (OUTFILE) || die "Closing $out: $!";
	print "Operator$opName extraction was succesful!\n";
}

sub smartDiff
{
	my ($opName, $raw, $oracle, $output) = @_;
	my @rowsRaw;
	my @rowsOracle;
	my @elemRaw;
	my @elemOracle;
	my %mapPos;
	
	open (FILE, "<$raw") || die "Opening $raw: $!";
	while(my $line = <FILE>) {
		next if($line !~ /^\d/);
		chomp($line);
		push @rowsRaw, $line;
	}
	close (FILE) || die "Closing $raw: $!";
	
	open (FILE, "<$oracle") || die "Opening $oracle: $!";
	while(my $line = <FILE>) {
		next if($line !~ /^\d/);
		chomp($line);
		push @rowsOracle, $line;
	}
	close (FILE) || die "Closing $oracle: $!";
	
	my @dimsRaw = split(' ', $rowsRaw[0]);
	my @dimsOracle = split(' ', $rowsOracle[0]);
	
	die "Unbalanced dimensions, Operator$opName matrix: $!" if($dimsRaw[0]  != $dimsOracle[0] || $dimsRaw[1] != $dimsOracle[1]);
	shift(@rowsRaw);
	shift(@rowsOracle);
	
	for(my $i = 0; $i < $dimsRaw[0]; $i++) {
		@elemRaw = split(' ', $rowsRaw[$i]);
		@elemOracle = split(' ', $rowsOracle[$i]);
		
		for(my $j = 0; $j < $dimsRaw[1]; $j ++) {
			if($elemRaw[$j] ne $elemOracle[$j]) {
				$mapPos{"($i,$j)"} = "$elemRaw[$j], $elemOracle[$j]";
			}
		}
	}
	
	if(scalar keys %mapPos) {
		open (FILE, ">$output") || die "Opening $output: $!";
		print FILE "Position    Raw    Oracle\n";
		print FILE "--------    ---    ------\n";
		foreach my $key (sort keys %mapPos) {
			print FILE "$key : $mapPos{$key}\n";
		}	
		close (FILE) || die "Closing $output: $!";
	}
		
	print "Smart diff for Operator$opName was successful.\n";
}

sub hookGprof
{
	my ($analysis, $arg) = @_;
	
	my $err = chdir($srcDir);
	die "Changing directory to $srcDir: $!" if(!$err);
	eval("system(\"gprof $arg\");");
	if($@) {
		my $subr = (caller(0))[3];
		die "$subr: $@";
	}
	$err = chdir($testDir);
	die "Changing directory to $testDir: $!" if(!$err);
	
	print "$analysis:Gprof was successful.\n" if($verbose);
	
}

sub hookDiff
{
	my ($analysis, $arg) = @_;
	
	eval("system(\"diff $arg\");");
	if($@) {
		my $subr = (caller(0))[3];
		die "$subr: $@";
	}

	print "$analysis:Diff was successful.\n" if($verbose);
	
}

sub hookGrep
{
	my ($analysis, $arg) = @_;
	
	eval("system(\"grep $arg\");");
	if($@) {
		my $subr = (caller(0))[3];
		die "$subr: $@";
	}
	
	print "$analysis:Grep was successful.\n" if($verbose);
}

sub moveFiles
{
	my @files = ("data$testNum.txt", "tst$testNum.txt", "SystemStackdata$testNum.txt", "EnvironStackdata$testNum.txt", "timeEvolution$testNum.txt");
	my $destination = $resultsDir;
	my $err = chdir($srcDir);
	die "Changing directory to $srcDir: $!" if(!$err);
	
	foreach my $f (@files) {
		if(validateFile($f)) {
			$err = system("mv $f $destination");
			die "Moving $f to $destination: $!" if($err);
		}
	}
	$err = chdir($testDir);
	die "Changing directory to $testDir: $!" if(!$err);
}

sub removeFiles
{
	my @files = ("Makefile*", "dmrg.*", "observe.*", "freeSystem*", "input.*", "raw$testNum.txt", "gmon.out");

	my $err = chdir($srcDir);
	die "Changing directory to $srcDir: $!" if(!$err);
	
	system("rm @files >& /dev/null");
	$err = chdir($testDir);
	die "Changing directory to $testDir: $!" if(!$err);
	system("rm @files >& /dev/null");
}
