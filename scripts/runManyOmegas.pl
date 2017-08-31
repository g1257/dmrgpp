#!/usr/bin/perl -w
#
use strict;

my ($n,$step)=@ARGV;

for (my $i=0;$i<$n;$i++) {
	my $omega = $i*$step;
	my $outputFile = runThisOmega($i,$omega);
	my $x = getResults($outputFile,$i,"OMEGA ".$omega);
	print "$omega $x\n";
}

sub runThisOmega
{
	my ($ind,$omega)=@_;
	my $inputFile = createInput($ind,$omega);
	my $outputFile = "output.out";
	system("./dmrg $inputFile >& $outputFile");
	return $outputFile;
}

sub getResults
{
	my ($f,$ind,$omega)=@_;
	open(FILE, "<", $f) or die "Cannot open $f: $!\n";
	my $ptrn = $omega;
	$ptrn =~ s/\./\\./;
	my $val=0;
	while(<FILE>) {
		if (/^$ptrn/) {
			my @temp = split;
			($temp[0]==$omega) or die "$temp[0]!=$omega\n";
			$val = $temp[1];
		}
	}
	close(FILE);
	return $val;
}

sub createInput
{
	my ($ind,$omega)=@_;
	my $inputFile = "input.inp";
	open(FOUT, ">", "$inputFile") or die "Cannot open $inputFile for writing: $!\n";
	print FOUT<<EOF;
TotalNumberOfSites=8 
NumberOfTerms=1
DegreesOfFreedom=1
GeometryKind=chain
GeometryOptions=ConstantValues
Connectors 1 1.0

hubbardU	8   0 0 0 0   0 0 0 0
potentialV	16  0 0 0 0   0 0 0 0 
		    0 0 0 0   0 0 0 0
SolverOptions=hasQuantumNumbers,wft,nosu2,DynamicTargetting
Version=6ce41a4b7dfa08978e53fa756f7f139e2fb18251
OutputFile=data.txt
InfiniteLoopKeptStates=200 
FiniteLoops 8   3 200 0 -3 200 0 -3 200 0 3 200 0 
                3 200 0 -3 200 0 -3 200 0 3 200 0

TargetQuantumNumbers 2 0.5 0.5

TSPFilename=tst.txt
TSPTau=0.1
TSPTimeSteps=5
TSPAdvanceEach=100
TSPSites 1 4 
TSPLoops 1 0 

TSPOperator=raw 
RAW_MATRIX 
4 4
0.0    0.0    0.0   0.0
1.0    0.0    0.0   0.0
0.0    0.0    0.0   0.0 
0.0    0.0    1.0   0.0 
FERMIONSIGN=-1
JMVALUES 0 0
AngularFactor=1
TSPOmega=$omega
TSPEta=0.05

EOF
	close(FOUT);
	return $inputFile;
}

