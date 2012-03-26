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

use strict;
use TestSuiteGlobals;

package TestSuiteDmrg;

#Custom routine that creates the dmrg executable, if necessary, and runs it
sub runDmrg
{
	my ($inputFile,$raw) = @_;

	#print STDERR "raw=$raw\n";
	die "Missing input file $inputFile: $!\n" unless (-r "$inputFile");
	
	my ($specFile, $specKey) = TestSuiteGlobals::getSpecFileAndKey();
	my $executable = $TestSuiteGlobals::srcDir."dmrg-".$specKey;

	createExecutable($specFile,$specKey,"dmrg") unless (-x "$executable");

	my $arg = "$executable $inputFile &> $raw";
# 	grep {s/&//} $arg if($verbose);
	
	print "Running dmrg test...\n";
	my $err = chdir($TestSuiteGlobals::srcDir);
	die "Changing directory to $TestSuiteGlobals::srcDir: $!" if(!$err);
	$err = system($arg);
	die "Running test using $executable with $inputFile: $!" if($err);
	$err = chdir($TestSuiteGlobals::testDir);
	die "Changing directory to $TestSuiteGlobals::testDir: $!" if(!$err);
	print "Completion of dmrg test.\n";
}

#Custom routine that creates the observe executable, if necessary, and runs it
sub runObserve
{
	my ($inputFile, $raw,$obsOptions) = @_;
	$obsOptions = "ccnnszsz" if (!defined($obsOptions));
	my ($specFile, $specKey) = TestSuiteGlobals::getSpecFileAndKey();
	my $executable = $TestSuiteGlobals::srcDir."observe-".$specKey;
	
	createExecutable($specFile,$specKey,"observe") unless (-x "$executable");
	
	my $arg = "$executable $inputFile $obsOptions &> $raw";
# 	grep {s/&//} $arg if($verbose);
	#die "******-> $arg\n";	
	print "Running observe test...\n";
	my $err = chdir($TestSuiteGlobals::srcDir);
	die "Changing directory to $TestSuiteGlobals::srcDir: $!" if(!$err);
	$err = system($arg);
	#die "Running test using $executable with $inputFile: $!" if($err);
	$err = chdir($TestSuiteGlobals::testDir);
	die "Changing directory to $TestSuiteGlobals::testDir: $!" if(!$err);
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
	
	my $err = chdir($TestSuiteGlobals::srcDir);
	die "Changing directory to $TestSuiteGlobals::srcDir: $!" if(!$err);
	print "Configuring $execType in Test $TestSuiteGlobals::testNum...\n";
	$err = system($arg1);
	die "Configuration error using $configFile with $specFile: $!" if($err);
	print "Creating $execType executable for Test $TestSuiteGlobals::testNum...\n";
	$err = system($arg2);
	die "Make command for $execType: $arg2\n" if($err);

	my $executable= $execType."-".$refKey;
	$err = rename($execType, $executable);
	die "Renaming $execType to $executable: $!" if(!$err);
	
	$err = chdir($TestSuiteGlobals::testDir);
	die "Changing directory to $TestSuiteGlobals::testDir: $!" if(!$err);
	
	print "\u$execType executable was succesfully created.\n";

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
	
	if($dimsRaw[0]  != $dimsOracle[0] || $dimsRaw[1] != $dimsOracle[1]) {
		print STDERR "$0: smartDiff: Unbalanced dimensions, Operator$opName matrix.";
		print STDERR " SmartDiff will be omitted\n";
		return;
	}

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


1;
