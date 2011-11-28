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
# use Cwd 'abs_path';
package TestSuiteGlobals;

my ($testDir, $srcDir,$inputsDir);

sub init
{
# my $PATH = $testDir = $srcDir = Cwd::abs_path($0);
# chomp(my $filename = `basename $0`);
# $testDir =~ s/$filename$//;
# $srcDir =~ s/TestSuite.*/src\//;
$TestSuiteGlobals::testDir=`pwd`;
chomp($TestSuiteGlobals::testDir);
$TestSuiteGlobals::testDir.="/";
$TestSuiteGlobals::srcDir=$TestSuiteGlobals::testDir."../src/";
$TestSuiteGlobals::inputsDir = $TestSuiteGlobals::testDir."inputs/";
print "$TestSuiteGlobals::inputsDir\n";
}

#Hook routine for the 'gprof' command
sub hookGprof
{
	my ($analysis, $arg) = @_;
	
	my $err = chdir($TestSuiteGlobals::srcDir);
	die "Changing directory to $TestSuiteGlobals::srcDir: $!" if(!$err);
	eval("system(\"gprof $arg\") == 0 || die;");
	if($@) {
		my $subr = (caller(0))[3];
		die "$subr: $@";
	}
	$err = chdir($TestSuiteGlobals::testDir);
	die "Changing directory to $TestSuiteGlobals::testDir: $!" if(!$err);
	
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

#Retrieves model spec file of current test and returns the file with its hash key
sub getSpecFileAndKey
{
	my $specFile = $TestSuiteGlobals::inputsDir."model$TestSuiteGlobals::testNum.spec";
	
	my $specKey = substr(`md5sum $specFile`,0,8);
	$specKey .= substr(`git rev-parse HEAD`,0,4);
	
	return $specFile, $specKey;
}
1;
