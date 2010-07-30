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

my ($testNumber, $lastTest) = ("","");
my ($all, $verbose, $rmFlag,$noModel, $observeFlag) = (0,0,0,0,0,0);
my ($testDir, $srcDir);
my $PATH = $testDir = $srcDir = abs_path($0);
chomp(my $filename = `basename $0`);
$testDir =~ s/$filename$//;
$srcDir =~ s/TestSuite.*/src\//;
my $oraclesDir = $testDir."oracles/";	
my $inputsDir = $testDir."inputs/";
my $globalRunFlag = 0;
my $executable = "";

$SIG{INT} = \&exit_handler;
$SIG{__WARN__} = sub {die @_};

sub exit_handler
{
	print "\nOracle creator aborted -> Manual cancellation\n";
	cleanUp();
	kill 9, $$;
}

sub cleanUp
{
	#Additional files can be added to @files to be removed
	my @files = ("freeSystem*", "input.*", "raw$testNumber.txt", "gmon.out", "data$testNumber.txt", "tst$testNumber.txt", "SystemStackdata$testNumber.txt", "EnvironStackdata$testNumber.txt", "timeEvolution$testNumber.txt", "stderrAndOut$testNumber.txt");

	my $err = chdir($srcDir);
	die "Changing directory to $srcDir: $!" if(!$err);
	
	system("rm @files >& /dev/null");
	$err = chdir($testDir);
	die "Changing directory to $testDir: $!" if(!$err);
	system("rm @files >& /dev/null");
}

eval {
	die $! if(!GetOptions("n=i" => \$testNumber, "l=i" => \$lastTest, "all" => \$all,
	"verbose" => \$verbose, "remove" => \$rmFlag, "manual" => \$noModel, "observe" => \$observeFlag));
	
	$globalRunFlag = 1;
	
	while($globalRunFlag) {
		runScript();
	}
};
if($@) {
	print "\nOracle creator aborted -> $@";
	cleanUp();
	kill 9, $$;
}

sub runScript
{
	print "*******INITIAL CONFIGURATIONS*******\n";
	die "$!" if(!validateDirectory($srcDir));
	die "$!" if(!validateDirectory($testDir));
	die "$!" if(!validateDirectory($oraclesDir));
	
	if(!($all)) {
		if($testNumber eq "") {
			selectTest();
		} else {
			if(!validateTest(getAvailableTests())) {
				selectTest();
			}
		}
	}			
	
	if($all) {
		runAllOracles(0);
	} elsif($testNumber < 0) {
		runAllOracles(-$testNumber);
	} else {
		runOracle();
	}
	
	print "*******FINAL CONFIGURATIONS*******\n";
	
	#Stop the testsuite program
	$globalRunFlag = 0;
}

#Displays available tests until user selects a valid one
sub selectTest
{
	my $available = getAvailableTests();
	
	while() {
		print "Type the number of the test you want to run.\n";
		print "Available tests: $available\n";
		print "Default is 0 (press ENTER): ";
		chomp($testNumber = <STDIN>);
		$testNumber = 0 if(!$testNumber);
		last if(validateTest($available));
	}
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

#Verifies if the test selected is valid
sub validateTest
{
	my ($available) = @_;
	
	if($testNumber =~ /(^-?\d+$)/) {
		my $searchNum = abs($testNumber);
		my @found = grep(/$searchNum/, split(/ /,$available));
		return 1 if(@found);
	}
	
	print "\n<Error>: An incorrect test was selected! Try again.\n\n";
	return 0;
}

#Verifies if a directory exists; creates a directory for the oracles if required
sub validateDirectory	
{
	my ($dir) = @_;
	
	if(-d $dir) {		
		return 1;
	} elsif($dir eq $oraclesDir){
		mkdir($dir) || die "Making directory $dir: $!";
		print "Directory created: $dir\n";
		return 1;
	}
	
	return 0;
}

#Verifies if a file exists; creates a hash table file if required
sub validateFile	
{
	my ($file) = @_;
	
	if(-e $file) {		
		return 1;
	} 
	
	return 0;
} 

sub runOracle
{
	my $configFile = "configure.pl";
	my $inputFile = $inputsDir."input$testNumber.inp";
	
	if(validateFile($inputFile)) {
		print "*******START OF TEST $testNumber*******\n";
		my $execType = "dmrg";
		$executable = createExecutables($configFile,$execType);
		runDmrg($inputFile);
		
		if($observeFlag)
		{
			$execType = "observe";
			$executable = createExecutables($configFile,$execType);
			runObserve($inputFile);
		}
		
		hookGprof();
		print "*******END OF TEST $testNumber*******\n";
	} else {
		die "$!";
	}
	
	print "Saving result files...\n";
	moveFiles();
	if($rmFlag) {
		print "Removing temporary files...\n";
		removeFiles() ;
	}	
	$executable = "";
	$testNumber = "";
}

sub runAllOracles
{
	my ($start) = @_;
	my @nonFunctionalTests = (24,41,42,60,104,105,106,107,108,109,110,111,124,125,141,142,160);
	my @testsList = split(/ /,getAvailableTests());
	
	if($lastTest ne "") {
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
		$testNumber = $testsList[$i];
		runOracle();
		last if($testsList[$i] == $lastTest);
	}
}

sub createExecutables
{
	my ($configFile,$execType) = @_;
	my $tempNum = $testNumber;
	$tempNum -= 100 if($testNumber >= 100);
	my $specFile = $inputsDir."model$tempNum.spec";
	my $arg1 = "./$configFile < $specFile >& /dev/null";
	my $arg2 = "make $execType -f Makefile >& /dev/null";

	grep {s/<.*//} $arg1 if($noModel);
	grep {s/>.*//} $arg1 if($verbose);
	grep {s/>.*//} $arg2 if($verbose);
	
	my $err = chdir($srcDir);
	die "Changing directory to $srcDir: $!" if(!$err);
	print "Configuring $execType in Test $testNumber...\n";
	$err = system($arg1);
	die "Configuration error using $configFile with $specFile: $!" if($err);
	print "Creating $execType executable for Test $testNumber...\n";
	$err = system($arg2);
	die "Make command for $execType: $!" if($err);
	$err = chdir($testDir);
	die "Changing directory to $testDir: $!" if(!$err);
	
	print "\u$execType executable was succesfully created.\n";
	
	return $srcDir.$execType;	
}

sub runDmrg
{
	my ($inputFile) = @_;
	my $raw = $oraclesDir."stderrAndOut$testNumber.txt";
	my $dataFile = $srcDir."data$testNumber.txt";
	my $arg = "$executable $inputFile >& $raw";
	
	die "Missing input file: $!" if(!validateFile($inputFile));
	grep {s/&//} $arg if($verbose);
	
	print "Running dmrg test...\n";
	my $err = chdir($srcDir);
	die "Changing directory to $srcDir: $!" if(!$err);
	$err = system($arg);
	die "Running test using $executable with $inputFile: $!" if($err);
	$err = chdir($testDir);
	die "Changing directory to $testDir: $!" if(!$err);
	print "Completion of dmrg test.\n";
	
	extractEnergy($dataFile);	
}

sub runObserve
{
	my ($inputFile) = @_;
	my $raw = $srcDir."raw$testNumber.txt";
	my $arg = "$executable $inputFile >& $raw";
	
	die "Missing input file: $!" if(!validateFile($inputFile));
	grep {s/&//} $arg if($verbose);
	
	print "Running observe test...\n";
	my $err = chdir($srcDir);
	die "Changing directory to $srcDir: $!" if(!$err);
	$err = system($arg);
	die "Running test using $executable with $inputFile: $!" if($err);
	$err = chdir($testDir);
	die "Changing directory to $testDir: $!" if(!$err);
	print "Completion of observe test.\n";
	
	extractOperatorC($raw);
	extractOperatorN($raw);
	extractOperatorSz($raw);
}

sub extractEnergy
{
	my ($data) = @_;
	my $energyOut = $oraclesDir."e$testNumber.txt";
        my $arg = "grep Energy $data > $energyOut";
	my $err = system($arg);
	die "Grep stopped: $!" if($err);
        print "Energy extraction was successful.\n" if($verbose);
}

sub extractOperatorC
{
	my ($raw) = @_;
	my $cOut = $oraclesDir."operatorC$testNumber.txt";
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
	print "OperatorC extraction was successful.\n" if($verbose);
}


sub extractOperatorN
{
	my ($raw) = @_;
	my $nOut = $oraclesDir."operatorN$testNumber.txt";
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
	print "OperatorN extraction was successful.\n" if($verbose);
}


sub extractOperatorSz
{
	my ($raw) = @_;
	my $sOut = $oraclesDir."operatorSz$testNumber.txt";
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
	print "OperatorSz extraction was successful.\n" if($verbose);
}

sub hookGprof
{
	my $profFile = $oraclesDir."prof$testNumber.txt";
	my $arg = "gprof $executable > $profFile";
	
	my $err = chdir($srcDir);
	die "Changing directory to $srcDir: $!" if(!$err);
	$err = system($arg);
	die "Gprof stopped: $!" if($err);
	$err = chdir($testDir);
	die "Changing directory to $testDir: $!" if(!$err);
	
	print "Profiling was successful.\n" if($verbose);
}

sub moveFiles
{
	#Additional files can be added to @files and $destination can be modified as wanted
	my @files = ("data$testNumber.txt", "tst$testNumber.txt", "SystemStackdata$testNumber.txt", "EnvironStackdata$testNumber.txt", "timeEvolution$testNumber.txt");
	my $destination = $oraclesDir;
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

#Cleans and removes unnecessary and temporary files during the final processing phase
sub removeFiles
{
	#Additional files can be added to @files to be removed
	my @files = ("Makefile*", "dmrg.*", "observe.*", "freeSystem*", "input.*", "raw$testNumber.txt", "gmon.out", "tst$testNumber.txt", "SystemStackdata$testNumber.txt", "EnvironStackdata$testNumber.txt", "timeEvolution$testNumber.txt", "stderrAndOut$testNumber.txt");

	my $err = chdir($srcDir);
	die "Changing directory to $srcDir: $!" if(!$err);
	
	system("rm @files >& /dev/null");
	$err = chdir($testDir);
	die "Changing directory to $testDir: $!" if(!$err);
	system("rm @files >& /dev/null");
}