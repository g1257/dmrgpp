#!/usr/bin/perl 
=pod
// BEGIN LICENSE BLOCK
Copyright © 2009 , UT-Battelle, LLC
All rights reserved

[DMRG++, Version 2.0.0]

*********************************************************
THE SOFTWARE IS SUPPLIED BY THE COPYRIGHT HOLDERS AND
CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED
WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A
PARTICULAR PURPOSE ARE DISCLAIMED. 

Please see full open source license included in file LICENSE.
*********************************************************


// END LICENSE BLOCK
=cut
use warnings;
use strict;

sub askQuestions
{
	if ($model=~/febasedsc/i) {
		$geometry="ladderfeas";
	} else {
		print "What geometry do you want to use?\n";
		print "Available: 1D, Ladder or LadderFeAs\n";
		print "Default is: 1D (press ENTER): ";
		$_=<STDIN>;
		chomp;
		if ($_ eq "" or $_ eq "\n") {
			$_="1d";
		}
		$geometry=$_;
	}
	
	if ($geometry=~/ladder$/i) {
		print "Enter the leg of ladder\n";
		print "Available: any even number\n";
		print "Default is: 2 (press ENTER): ";
		$_=<STDIN>;
		chomp;
		if ($_ eq "" or $_ eq "\n") {
			$_=2;
		}
		$legOfLadder=$_;
	}
	
	print "What targetting do you want?\n";
	print "Available: GroundStateTargetting TimeStepTargetting\n";
	print "Default is: GroundStateTargetting (press ENTER): ";
	$_=<STDIN>;
	chomp;
	if ($_ eq "" or $_ eq "\n") {
		$_="GroundStateTargetting";
	}
	$targetting = $_;
	
	print "Enter the number of kept states for the infinite loop\n";
	print "Available: any\n";
	print "Default is: 64 (press ENTER): ";
	$_=<STDIN>;
	chomp;
	if ($_ eq "" or $_ eq "\n") {
		$_=64;
	}
	$infiniteKeptStates=$_;

	print "Enter the total number of sites (system+environment)\n";
	print "Available: any\n";
	print "Default is: 16 (press ENTER): ";
	$_=<STDIN>;
	chomp;
	if ($_ eq "" or $_ eq "\n") {
		$_=16;
	}
	$linSize=$_;

	askAboutFiniteLoops();		
	
	print "Enter the value of the $connectors\n";
	print "Available: any\n";
	print "Default is: 1.0 (press ENTER): ";
	$_=<STDIN>;
	chomp;
	if ($_ eq "" or $_ eq "\n") {
		$_="1.0";
	}
	$connectorValue=$_;
	
	if (defined($connectors2)) {
		print "Enter the value of the $connectors2\n";
		print "Available: any\n";
		print "Default is: 1.0 (press ENTER): ";
		$_=<STDIN>;
		chomp;
		if ($_ eq "" or $_ eq "\n") {
			$_="1.0";
		}
		$connectorValue2=$_;
	}
	
	if ($model=~/hubbard/i) {
		askQuestionsHubbard();
	} elsif ($model=~/febasedsc/i) {
		askQuestionsFeAs();
	} elsif ($model=~/tjoneorbital/i) {
#askQuestionsTjOneOrbital(); # FIXME: ask questions about potential
	}

	print "Do you want to run with SU(2) symmetry enabled?\n";
	print "Available: yes or no\n";
	print "Default is: no (press ENTER): ";
	$_=<STDIN>;
	chomp;
	if ($_ eq "" or $_ eq "\n") {
		$_="no";
	}
	$su2Symmetry=$_;
	if ($su2Symmetry=~/y/i) {
		print "There will be ".($linSize)." sites in total, so...\n";
		print "how many total electrons should I consider?\n";
		print "Default is: ".($linSize)." (press ENTER): ";
		$_=<STDIN>;
		chomp;
		if ($_ eq "" or $_ eq "\n") {
                	$_=$linSize;
        	}
		$electrons=$_;
		print "If there are ".($linSize)." sites in total...\n";
                print "What value of angular momentum j do you want?\n";
                print "Default is: 0 (press ENTER): ";
                $_=<STDIN>;
                chomp;
                if ($_ eq "" or $_ eq "\n") {
                        $_=0;
                }
                $momentumJ=$_;
	}

}


sub askQuestionsHubbard()
{ 	
	print "Enter the value of the Hubbard U values\n";
	print "Available: any\n";
	print "Default is: 1.0 (press ENTER): ";
	$_=<STDIN>;
	chomp;
	if ($_ eq "" or $_ eq "\n") {
		$_="1.0";
	}
	$hubbardUvalue=$_;
	
	print "Enter the value of the potential values\n";
	print "Available: any\n";
	print "Default is: 0.0 (press ENTER): ";
	$_=<STDIN>;
	chomp;
	if ($_ eq "" or $_ eq "\n") {
		$_="0.0";
	}
	$potentialVvalue=$_;
}

sub askQuestionsFeAs()
{
	print STDERR "Please check your input file is written correctly\n";
}

sub createInput
{
	system("cp input.inp input.bak") if (-r "input.inp");
	open(FOUT,">input.inp") or die "Cannot open file input.inp for writing: $!\n";
	my $connectorValues=createConnectors($connectorValue);
	my $version=getVersion();
	my $qns = "2 0.5 0.5";
	my $inputForTimeEvolution = getTimeEvolutionInput();
	my $geometryName = getGeometryName();
	
	if ($su2Symmetry=~/y/i) {
		$_ = $electrons/(2*$linSize);
		$momentumJ /= ($linSize);
		$qns="3 $_ $_ $momentumJ\n";
		$su2Symmetry="useSu2Symmetry";
	} else {
		$su2Symmetry="nosu2";
	}
	my $terms = 1; # it is 2 for t-j model
	my $edof = 1;
	print FOUT "TotalNumberOfSites=$linSize\n";
	print FOUT "NumberOfTerms=$terms\n";
	print FOUT "DegreesOfFreedom=$edof\n";
	print FOUT "GeometryKind=$geometryName\n";
	print FOUT "GeometryOptions=ConstantValues\n";
	print FOUT "Connectors 1 1.0\n"; # FIXME only valid for hubbard model
	if ($model=~/febasedsc/i) {
	$qns = "3 1.0 1.0 0.0"; 
	print FOUT<<EOF;
hoppings 2 2
-0.058 0
0 -0.2196
hoppings 2 2
-0.2196 0
0 -0.058
hoppings 2 2
+0.20828 +0.079
+0.079 +0.20828
hoppings 2 2
+0.20828 -0.079
-0.079 +0.20828
EOF
	} else {
		# FIXME
	}
	
	if ($model=~/hubbard/i or $model=~/febasedsc/i) {	
		my $hubbardu=createHubbardU();
		print FOUT<<EOF;
hubbardU	$hubbardu
EOF
	}
	if ($model=~/hubbard/i or $model=~/febasedsc/i or $model=~/tjoneorbital/i) {
		my $potentialv=createPotentialV();

print FOUT<<EOF;
potentialV	$potentialv
density=1.0
EOF
	}

	my $hasThreads="";
	$hasThreads = "hasThreads" if ($pthreads);

	print FOUT<<EOF;
SolverOptions=hasQuantumNumbers,nowft,$su2Symmetry,$hasLoops,$hasThreads,$targetting
Version=$version
OutputFile=data.txt
InfiniteLoopKeptStates=$infiniteKeptStates
FiniteLoops $finiteLoops
TargetQuantumNumbers $qns
   
EOF
	print FOUT "Threads=$nthreads\n" if ($pthreads);
	print FOUT "$inputForTimeEvolution\n\n" if ($targetting=~/timestep/i);
	print STDERR "File input.inp has been written\n";
	close(FOUT);
}

sub getTimeEvolutionInput
{
	return "#NO_TIME_EVOLUTION" unless ($targetting=~/timestep/i);
	my $ret = <<EOF;
FILENAME tst.txt
TIMESTEP 0.1
MAXTIMES 4 
ADVANCEEACH 5 
SITE  2 10 11 
STARTINGLOOP 2 0 0 

TIMEEVOLUTION raw 
RAW_MATRIX 4 4
0.0    0.0    0.0   0.0
0.0   0.0    0.0   -1.0
0.0    0.0    0.0   1.0 
0.0    0.0    0.0   0.0 
FERMIONSIGN -1
JMVALUES 0 0
angularFactor 1

TIMEEVOLUTION raw
RAW_MATRIX 4 4
0.0    0.0    0.0   0.0
1.0    0.0    0.0   0.0
1.0    0.0    0.0   0.0
0.0    0.0    0.0   0.0
FERMIONSIGN -1
JMVALUES 0 0
angularFactor 1
EOF
	return $ret;
}

sub createConnectors
{
	my ($value)=@_;
	if ($geometry=~/ladder$/i) {
		return createConnectorsLadder($value);
	} elsif ($geometry=~/1d/i) {
		return createConnectorsChain($value);
	}
	return "UNDEFINED_GEOMETRY"; # should never reach here
}

sub createConnectorsChain
{
	my ($value)=@_;
	my $sbOpen=""; #"["
	my $sbClose=""; #"]"
	my $spacer=" "; #",";
	my $ret=$sbOpen;
	my ($i,$j);
	for ($j=0;$j<$linSize;$j++) {
		$ret=$ret.$sbOpen;
		for ($i=0;$i<$linSize;$i++) {
			if (abs(int($i-$j))==1) {
				$ret=$ret.$value;
			} else {
				$ret=$ret."0.0";
			}
			$ret=$ret.$spacer unless ($i==$linSize-1);
		}
		$ret=$ret."$sbClose\n";
		$ret=$ret.$spacer if ($j<$linSize-1);
	}
	$ret=$ret.$sbClose unless ($sbClose eq "");
		
	return $ret;
}

sub createConnectorsLadder
{
	my ($value)=@_;
	my ($x,$y,$i,$j);
	my @matrix;
	for ($x=0;$x<$linSize;$x++) {
		for ($y=0;$y<$linSize;$y++) {
			$matrix[$x][$y]=0;
		}
	}
	
	for ($x=0;$x<int($linSize/$legOfLadder);$x++) {
		for ($y=0;$y<$legOfLadder;$y++) {
			$i = $y +$x*$legOfLadder;
			if ($y+1<$legOfLadder) {
				$j=$y+1 + $x*$legOfLadder;
				$matrix[$i][$j]=$matrix[$j][$i]=$value;
			}
			if ($x+1<int($linSize/$legOfLadder)) {
				$j=$y + ($x+1)*$legOfLadder;
				$matrix[$i][$j]=$matrix[$j][$i]=$value;
			}
		}
	}
	my $sbOpen=""; #"["
	my $sbClose=""; #"]"
	my $spacer=" "; #",";
	my $ret=$sbOpen;
	for ($x=0;$x<$linSize;$x++) {
		$ret=$ret.$sbOpen;
		for ($y=0;$y<$linSize;$y++) {
			$ret=$ret.$matrix[$x][$y];
			$ret=$ret.$spacer unless ($y==$linSize-1);
		}
		$ret=$ret."$sbClose\n";
		$ret=$ret.$spacer unless ($x==$linSize-1);
	}
	$ret=$ret.$sbClose unless ($sbClose eq "");
	
	return $ret;
}	
	
sub createHubbardU
{
	#my $ret="[";
	return "4 0.0 0.0 0.0 0.0" if ($model=~/febasedsc/i);
	my $ret="$linSize ";
	my ($j);
	for ($j=0;$j<$linSize-1;$j++) {
		$ret=$ret."$hubbardUvalue ";
	}
	$ret=$ret."$hubbardUvalue"; #]";
	return $ret;
}

sub createPotentialV
{
	#my $ret="[";
	return "4 0.0 0.0 0.0 0.0" if ($model=~/febasedsc/i);
	my $ret=" ".(2*$linSize)." ";
	my ($j);
	for ($j=0;$j<2*$linSize-1;$j++) {
		$ret=$ret."$potentialVvalue ";
	}
	$ret=$ret."$potentialVvalue"; #]";
	return $ret;
}

sub getVersion
{
	my $version="NOGIT";
	my $hasGit=0;
	# First try to see if git is available
	my $tmp=system("git log -1 >& /dev/null");
	if ($tmp==0) {
		$hasGit=1;
	}
	$version=getGitVersion() if ($hasGit);
	
	 
	#If not then read it from a file
	if ($version=~/NOGIT/i) {
		if (open(VERSIONF,"version.txt")) {
			while(<VERSIONF>) {
				chomp;
				if (/commit (.*$)/) {
					$version=$1;
				}
			}
			close(VERSIONF);
		} else {
			$version="UNDEFINED";
			print STDERR "$0: (WARNING): Could not open version.txt for reading: $!\n";
		}
	}
	return $version;
}


sub getGitVersion
{
	my $version="NOGIT";
	open(PIPE,"git log -1 |") or return $version;
	while(<PIPE>) {
		chomp;
		if (/commit (.*$)/) {
			$version=$1;
		}
	}
	close(PIPE);
	if (open(VERSIONF,">version.txt")) {
		print VERSIONF "commit $version\n";
		close(VERSIONF);
	} else {
		print STDERR "$0: (WARNING): Could not open version.txt for writing: $!\n";
	}
	return $version;

}

sub askAboutFiniteLoops
{
	print "Do you want to do finite loops?\n";
        print "Available: y or n\n";
        print "Default is: y (press ENTER): ";
        $_=<STDIN>;
        chomp;
        if ($_ eq "" or $_ eq "\n") {
                $_="y";
        }
	$finiteLoops = "";
	$hasLoops="";
	if ($_=~/n/i) {
		$finiteLoops="1 1 100 0";
		$hasLoops="nofiniteloops";
		return;
	}
	my ($log,$m);
	my $x = $linSize/2-1;
	my $counter=1;
	while (1) {
		$m = addFiniteLoop();
		print "Do you want to add another finite loop?\n";
        	print "Available: y or n\n";
        	print "Default is: n (press ENTER): ";
        	$_=<STDIN>;
        	chomp;
        	if ($_ eq "" or $_ eq "\n") {
                	$_="n";
        	}
		$log = 0;
		last if ($_=~/n/i);
		$log =0;
		$finiteLoops = $finiteLoops." $x $m $log -$x $m $log -$x $m $log $x $m $log ";
		$counter++;
	}
	$finiteLoops = (4*$counter)." ".$finiteLoops." $x $m $log -$x $m $log -$x $m $log $x $m $log ";
}

sub addFiniteLoop
{
	print "What's the m for this loop?\n";
	print "Available: any\n";
	print "Default is: 60 (press ENTER): ";
	$_=<STDIN>;
	chomp;
	if ($_ eq "" or $_ eq "\n") {
		$_=60;
	}
	return $_;

}


sub getGeometryName
{
	my $geometryName = "UNKNOWN";
	if ($geometry=~/1d/i) {
		$geometryName = "chain";
	} elsif ($geometry=~/ladder$/i) {
		$geometryName = "ladder";
	} elsif ($geometry=~/ladderfeas/i) {
		$geometryName = "ladderx";
	}
	return $geometryName;
}
