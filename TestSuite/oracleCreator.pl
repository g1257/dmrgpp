#! /usr/bin/perl -w

use strict;
use Getopt::Long;
my ($testNumber, $all, $verbose, $rmFlag,$noModel, $observeFlag) = (undef,0,0,0,0,0);

#User can either select a specific test number, a range of tests (negative number), or all tests
#When clean is enabled all temporary files created will be deleted
GetOptions("n=i" => \$testNumber, "all" => \$all, "v" => \$verbose, "rm" => \$rmFlag
, "m" => \$noModel, "o" => \$observeFlag);

my $oraclesDir = "oracles/";		#Specifies the output folder for the oracles results
my $srcDir = "../src/";
my $testDir = "../TestSuite/";

selectTest();

sub selectTest
{
	if(validateDirectory($srcDir) && validateDirectory($testDir) && validateDirectory($oraclesDir)) {
		if(!defined($testNumber) && !($all)) {
			$testNumber = selectOracle();
		}
		
		
		if($all) {
			runAllOracles(0);
		} elsif($testNumber < 0) {
			runAllOracles(-$testNumber);
		} else {
			runOracle($testNumber);
		}
	}
}

#Verify if output directory exists in the hierarchy
#If not it is created
sub validateDirectory	
{
	my ($dir) = @_;
	
	if(!-d $dir) {
		if($dir eq $oraclesDir) {
			my $rCode = system("mkdir $dir");
			die "Error making $dir: $!\n" if($rCode);
			print "Directory created: $dir\n";
		} else {
			die "Error searching for $dir: $!\n";
		}
	}
	
	return 1;
}

sub selectOracle
{
	my $available = getAvailableTests();	#Get a numerical list of the available tests
	my @temp;
	my @testArray;
	my $searchNum;
	#Here the user selects a test to be run
	while() {
		print "Type the test number for oracles creation.\n";
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
	my $descriptionFile = "inputs/descriptions.txt";
	
	open(FILE,$descriptionFile) || die "Error opening $descriptionFile: $!\n";
	while(<FILE>) {
		last if(/^\#TAGEND/);		#TAGEND delimits the available tests in this document
		if(/(^\d+)\)/) {
			$available = $available.$1." ";	
		}
		print if($verbose);	#Display the available tests with their descriptions provided
	}
	close(FILE);
	print "\n" if($verbose);

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
	} else {
		print "Error: $file does not exists.\n";
		return 0;
	}
}

sub runOracle
{
	my ($tn) = @_;
	my $inputFile = $testDir."inputs/input$tn.inp";
	my $executable = $srcDir."dmrg";
	my $temp = $tn;
	$temp -= 100 if($tn >= 100);
	my $specFile = $testDir."inputs/model$temp.spec";
	
	exit() if(!validateFile($inputFile) || !validateFile($specFile));
	print "Preparing to run oracle of test $tn...\n";
	
	chdir("$srcDir");
	my $profFile = $testDir.$oraclesDir."prof$tn.txt";
	my $rCode = 0;
	
	if($noModel) {
		$rCode = system("./configure.pl");
		die "Error with configure.pl: $!\n" if($rCode);
	} else {
		$rCode = system("./configure.pl < $specFile");
		die "Error with configure.pl: $!\n" if($rCode);
	}
	
	$rCode = system("make -f Makefile");
	die "Error with the Makefile: $!\n" if($rCode);
	$rCode = system("./$executable $inputFile");
	die "Error running test: $!\n" if($rCode);
	profile($executable,$profFile);
	chdir("$testDir");
	
	my $dataFile = $srcDir."data$tn.txt";
	my $rawFile = "raw$tn.txt";
	my $energyOut = "oracles/e$tn.txt";
	my $tstFile = $srcDir."tst$tn.txt";
	my $envStack = $srcDir."EnvironStackdata$tn.txt";
	my $sysStack = $srcDir."SystemStackdata$tn.txt";
	my $line;
	
	open FILE,"<$specFile" || die "Error opening file: $!\n";
	while($line = <FILE>) {
		if(grep (/^TimeStepTargetting/, $line)) {
			$rCode = system("mv $tstFile oracles/");
			die "Error moving time step targetting: $!\n" if($rCode);
		} elsif(grep (/^DiskStack/,$line)) {
			$rCode = system("mv $envStack $sysStack oracles/");
			die "Error moving stack file: $!\n" if($rCode);
		}
	}
	close FILE;
	
	extractEnergy($dataFile,$energyOut);
	observables($tn,$inputFile,$dataFile,$rawFile) if($observeFlag);	
	$rCode = system("mv $dataFile oracles/");
	die "Error moving $dataFile: $!\n" if($rCode);
	removeFiles($tn) if($rmFlag);
	print "Completed running oracle for test $tn.\n";
	
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
	close (OUTFILE);	
	print "OperatorSz extraction was succesful!\n";
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

	my $cOut = $oraclesDir."operatorC$tn.txt";
	my $nOut = $oraclesDir."operatorN$tn.txt";
	my $sOut = $oraclesDir."operatorS$tn.txt";
	
	extractOperatorC($raw,$cOut);
	extractOperatorN($raw,$nOut);
	extractOperatorS($raw,$sOut);
	
}

sub profile
{
	my ($exe,$prof) = @_;
	my $rCode = system("gprof $exe > $prof"); 
	die "Error with gprof: $!\n" if($rCode); 
	 print "Profiling was succesful!\n";
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

sub runAllOracles
{
	my ($start) = @_;
	print "Preparing to run all oracles starting from test $start..\n";
	
	my $available = getAvailableTests();
	my @temp = split(/ /,$available);
	
	for (my $i=0;$i<$#temp+1;$i++) {
		next if ($temp[$i] eq "");
		next if ($temp[$i]<$start);
		runOracle($temp[$i]);
	}
	
	print "All oracles were done succesfully.\n";

}


sub extractEnergy
{
	my ($data,$energyOut) = @_;
        my $rCode = system("grep Energy $data > $energyOut");
	die "Error extracting energy values: $!\n" if($rCode);
        print "Energy extraction was succesful!\n";
}