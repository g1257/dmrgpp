#! /usr/bin/perl -w

use strict;
use Getopt::Long;
my ($testNumber, $all);
my $rmFLAG = 0;
GetOptions("n=i" => \$testNumber, "all" => \$all);

if(!defined($testNumber) && !defined($all)) {
	$testNumber = selectOracle();
}

if(defined($all)) {
	runAllOracles(0);
} elsif($testNumber < 0) {
	runAllOracles(-$testNumber);
} else {
	runOracle($testNumber);
}

sub selectOracle
{
	my $available = getAvailableTests(1);	#Get a numerical list of the available tests
	my @temp;
	my @testArray;

	#Here the user selects a test to be run
	while() {
		print "Type the test number for oracles creation.\n";
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
			print "\nERROR: An incorrect test was selected! Try again.\n\n";
		} else {
			print "Preparing to run oracle of test number $tn\n";
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

sub runOracle
{
	my ($testNum) = @_;
	
	chdir("../src/");
	system("perl configure.pl < ../TestSuite/inputs/model$testNum.spec"); # die "Error configure\n";
	system("make -f Makefile"); #die "Error make\n";
	system("./dmrg ../TestSuite/inputs/input$testNum.inp"); # die "Error dmrg\n";
	chdir("../TestSuite/");
	
	my $exe = "../src/dmrg";
	my $dataFile = "data$testNum.txt";
	my $rawOut = "raw$testNum.txt";
	my $inputFile = "inputs/input$testNum.inp";
	profile($exe,$testNum);
	cookRaw($testNum,$inputFile,$dataFile,$rawOut);
	extractEnergy($testNum,$dataFile);
	extractOperatorC($testNum,$rawOut);
	extractOperatorN($testNum,$rawOut);
	extractOperatorS($testNum,$rawOut);
	
	removeFiles($testNum) if($rmFLAG);
	print "Completed running oracle for test $testNum.\n";
	
}

sub extractOperatorC
{
	my ($tn,$file) = @_;
	my $output = "oracles/OperatorC$tn.txt";
	
	
}


sub extractOperatorN
{
	my ($tn,$file) = @_;
	my $output = "oracles/OperatorN$tn.txt";
}


sub extractOperatorS
{
	my ($tn,$file) = @_;
	my $output = "oracles/OperatorS$tn.txt";
}

sub cookRaw
{
	my ($tn,$input,$data,$output) = @_;
	
	chdir("../src/");
	system("make observe -f Makefile"); #die "Error\n";
	chdir("../TestSuite/");
	
	system("../src/observe $input $data > $output"); #die "Error\n";
	
}

sub profile
{
	my ($exe,$testNum) = @_;
	
	system("gprof $exe > prof$testNum.txt"); # die "Error\n";
	#cook the output
}

sub removeFiles
{
	#delete files with $rmFLAG
	
	
}

sub runAllOracles
{
	my ($start) = @_;
	print "Run all oracles starting from test #$start\n";
	
	my $available = getAvailableTests(0);
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
	my ($fnum,$fileAnalyzed) = @_;

        system("grep Energy $fileAnalyzed > oracles/etest.txt");

        print "Energy extraction from $fileAnalyzed is complete.\n";
}
