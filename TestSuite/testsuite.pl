#! /usr/bin/perl

=pod
// BEGIN LICENSE BLOCK
/*
Copyright (c) 2008-2011, UT-Battelle, LLC
All rights reserved

[DMRG++, Version 2.0.0]
[by G.A., Oak Ridge National Laboratory]
[TestSuite by E.P., Puerto Rico and ORNL]

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
my ($testDir, $srcDir);
my $PATH = $testDir = $srcDir = abs_path($0);
chomp(my $filename = `basename $0`);
$testDir =~ s/$filename$//;
$srcDir =~ s/TestSuite.*/src\//;
my $oraclesDir = $testDir."oracles/";	
my $resultsDir = $testDir."results/";
my $inputsDir = $testDir."inputs/";
my %dmrgHash = ();
my %observeHash = ();


#Get command line options
my $all = 0;
my $testNum;
my $lastTest;
die $! if(!GetOptions("n=i" => \$testNum, "l=i" => \$lastTest, "all" => \$all));

$testNum = selectTest() if (!defined($testNum));

if (!$all) {
	testSuite($testNum);
	exit(0);
}

if($testNum < 0) {
	$testNum = -$testNum;
}

runAllTests($testNum,$lastTest);

#Displays available tests until user selects a valid one
sub selectTest
{
	my $available = getAvailableTests();

	print "Type the number of the test you want to run.\n";
	print "Available tests: $available\n";
	print "Default is 0 (press ENTER): ";
	my $testNum = <STDIN>;
	chomp($testNum);
	return $testNum;
}

#Extracts from the descriptions file all the available tests that can be run
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
#		print if($verbose);	
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

#Runs a single test
sub testSuite
{
	my ($testNum)=@_;
	#$tempNum -= 100 if($testNum >= 100);
	my $procFile = $inputsDir."processing$testNum.txt";
	my $procLib = $inputsDir."processingLibrary.txt";
	
	if(-r $procLib) {
		if(-r "$procFile") {
			print "*******START OF TEST $testNum*******\n";
			my @analyses = extractAnalyses($procFile) ;
			(@analyses) ? (processing(@analyses, $procLib)) : (print "Test $testNum does not includes any processing analyses.\n");
			print "*******END OF TEST ".$testNum."*******\n";
		} else {
			die "Could not find $procFile: $!";
		}
	} else {
		die "$!";
	}
}

#Runs multiple tests by invoking the testSuite routine for each test
sub runAllTests
{
	my ($start,$lastTest) = @_;
	#All test numbers declared in @nonFunctionalTests will be skipped and not ran
	#test102 features wft+su2 which is currently broken
	my @nonFunctionalTests = (41,42,60,102,104,105,106,107,108,109,110,111,124,125,141,142,160);
	my @testsList = split(/ /,getAvailableTests());
	
	if(defined($lastTest)) {
		die "<Error>: Invalid tests range [$start,$lastTest].\n" if($lastTest < $start);
		print "Preparing to run all tests from Test $start to Test $lastTest.\n";
	} else {
		print "Preparing to run all tests starting from Test $start...\n";
		$lastTest = $testsList[$#testsList];
	}
	
	for (my $i=0;$i<=$#testsList;$i++) {
		next if ($testsList[$i] eq "");
		next if ($testsList[$i]<$start);
		next if(grep {$_ eq $testsList[$i]}@nonFunctionalTests);
		my $testNum = $testsList[$i];
		testSuite($testNum);
		last if($testsList[$i] == $lastTest);
	}
}

#Reads all analyses described in the processing file for a specific test
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

#Extracts the commands for each analysis of the current test, sends them to a key:value parser routine,
#and then the commands are sent to an interpreter routine for execution
sub processing
{
	my $lib = pop(@_);
	my (@analyses) = @_;
	my %procHash;
	my @commands;
	
	print "Starting processing phase...\n";
	
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
	
	#The following section of code resolves dependencies among the different analysis in the
	#processing file for the current test
	my %procCount;
	my %tmpHash = %procHash;
	my @dependencies;
	my $depFlag;
	my @depKeys;
	my $keyword = "CallOnce";	#Keyword that establishes dependencies in the processing library
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
	
	print "Ending processing phase...\n";
}

#Resolves the variables in the commands for each analysis of the current test
sub keyValueParser
{
	my ($opsRef) = @_;
	my @tmpKV;
	my $keyword = "Let";	#Keyword that marks key:value pairs used in the parsing action
	my %tmpHash;
	my @commands;
	my %varHash;
	#Only variables found in @varArray can be used in the processing library for substitution purposes
	my @varArray = ("\$executable", "\$srcDir", "\$inputsDir", "\$resultsDir", "\$oraclesDir", "\$testNum");
	#Variables in @nonSubsArray will not be resolved during the parsing, instead
	#they are resolved prior to executing the command
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

#Executes commands from the processing library by hooking the command to its routine
sub commandsInterpreter
{
	my $analysis = pop(@_);
	my (@commands) = @_;
	#Keywords in @metaLang are used to hook a command with its routine
	#Routines, in conjunction with meta language keywords, can be added to expand the runable commands in the processing library
	#The commands will be executed following the order below from left to right
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

#Subroutine for the 'Execute' keyword
#'Execute' is used when the user defines custom routines 
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
	
#	print "[$analysis]:Execute command was successful.\n" if($verbose);
}

#Custom routine that creates the dmrg executable, if necessary, and runs it
sub runDmrg
{
	my ($inputFile,$raw) = @_;

	die "Missing input file: $!" unless (-r "$inputFile");
	
	my ($specFile, $specKey) = getSpecFileAndKey();
	my $executable = $srcDir."dmrg-".$specKey;

	createExecutable($specFile,$specKey,"dmrg") unless (-x "$executable");

	my $arg = "$executable $inputFile &> $raw";
# 	grep {s/&//} $arg if($verbose);
	
	print "Running dmrg test...\n";
	my $err = chdir($srcDir);
	die "Changing directory to $srcDir: $!" if(!$err);
	$err = system($arg);
	die "Running test using $executable with $inputFile: $!" if($err);
	$err = chdir($testDir);
	die "Changing directory to $testDir: $!" if(!$err);
	print "Completion of dmrg test.\n";
}

#Custom routine that creates the observe executable, if necessary, and runs it
sub runObserve
{
	my ($inputFile, $raw,$obsOptions) = @_;
	$obsOptions = "ccnnszsz" if (!defined($obsOptions));
	my ($specFile, $specKey) = getSpecFileAndKey();
	my $executable = $srcDir."observe-".$specKey;
	
	createExecutable($specFile,$specKey,"observe") unless (-x "$executable");
	
	my $arg = "$executable $inputFile $obsOptions &> $raw";
# 	grep {s/&//} $arg if($verbose);
	#die "******-> $arg\n";	
	print "Running observe test...\n";
	my $err = chdir($srcDir);
	die "Changing directory to $srcDir: $!" if(!$err);
	$err = system($arg);
	die "Running test using $executable with $inputFile: $!" if($err);
	$err = chdir($testDir);
	die "Changing directory to $testDir: $!" if(!$err);
	print "Completion of observe test.\n";
}

#Configures the current test, either manually or automatically (model spec file), and creates the executable
sub createExecutable
{
	my ($specFile,$refKey, $execType) = @_;
	my $configFile = "configure.pl";
	my $arg1 = "./$configFile < $specFile &> /dev/null";
	my $arg2 = "make $execType -f Makefile &> /dev/null";

# 	grep {s/>.*//} $arg1 if($verbose);
# 	grep {s/>.*//} $arg2 if($verbose);
	
	my $err = chdir($srcDir);
	die "Changing directory to $srcDir: $!" if(!$err);
	print "Configuring $execType in Test $testNum...\n";
	$err = system($arg1);
	die "Configuration error using $configFile with $specFile: $!" if($err);
	print "Creating $execType executable for Test $testNum...\n";
	$err = system($arg2);
	die "Make command for $execType: $!" if($err);

	my $executable= $execType."-".$refKey;
	$err = rename($execType, $executable);
	die "Renaming $execType to $executable: $!" if(!$err);
	
	$err = chdir($testDir);
	die "Changing directory to $testDir: $!" if(!$err);
	
	print "\u$execType executable was succesfully created.\n";

}

#Retrieves model spec file of current test and returns the file with its hash key
sub getSpecFileAndKey
{
	my $specFile = $inputsDir."model$testNum.spec";
	
	my $specKey = substr(`md5sum $specFile`,0,8);
	$specKey .= substr(`git rev-parse HEAD`,0,4);
	
	return $specFile, $specKey;
}

#Retrieves the data for the Operators C, N, and Sz from the observe run
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
	#die "Here <----------- $opName $raw\n";
	print OUTFILE $op;
	close (OUTFILE) || die "Closing $out: $!";
#	print "Operator$opName extraction was successful.\n" if($verbose);
}

#Searches for differences between the data in the operators oracles with the recently computed operators
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
	
	open (FILE, ">$output") || die "Opening $output: $!";
	if(scalar keys %mapPos) {
		print FILE "(Row,Col)   Raw    Oracle\n";
		print FILE "--------    ---    ------\n";
		foreach my $key (sort keys %mapPos) {
			print FILE "$key : $mapPos{$key}\n";
		}
	}
	close (FILE) || die "Closing $output: $!";
		
#	print "Smart diff for Operator$opName was successful.\n" if($verbose);
}

#Hook routine for the 'gprof' command
sub hookGprof
{
	my ($analysis, $arg) = @_;
	
	my $err = chdir($srcDir);
	die "Changing directory to $srcDir: $!" if(!$err);
	eval("system(\"gprof $arg\") == 0 || die;");
	if($@) {
		my $subr = (caller(0))[3];
		die "$subr: $@";
	}
	$err = chdir($testDir);
	die "Changing directory to $testDir: $!" if(!$err);
	
#	print "[$analysis]:Gprof command was successful.\n" if($verbose);
}

#Hook routine for the 'diff' command
sub hookDiff
{
	my ($analysis, $arg) = @_;
	
	eval("system(\"diff $arg\");");
	if($@) {
		my $subr = (caller(0))[3];
		die "$subr: $@";
	}

#	print "[$analysis]:Diff command was successful.\n" if($verbose);
}

#Hook routine for the 'grep' command
sub hookGrep
{
	my ($analysis, $arg) = @_;
	
	eval("system(\"grep $arg\") == 0 || die;");
	if($@) {
		my $subr = (caller(0))[3];
		die "$subr: $@";
	}
	
#	print "[$analysis]:Grep command was successful.\n" if($verbose);
}
