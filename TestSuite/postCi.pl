#!/usr/bin/perl

use strict;
use warnings;
use Getopt::Long qw(:config no_ignore_case);
use Ci;

my ($memory,$failed,$noSu2,$help,$workdir,$golddir,$ranges);
GetOptions(
'n=s' => \$ranges,
'memory=i' => \$memory,
'f' => \$failed,
'nosu2' => \$noSu2,
'h' => \$help,
'g=s' => \$golddir,
'w=s' => \$workdir);

if (defined($help)) {
	print "USAGE: $0 [options]\n";
	print "\tIf no option is given examines all runs\n";
	print Ci::helpFor("-n");
	print Ci::helpFor("-w");
	print "\t-g golddir\n";
	print "\t\tUse golddir for oracles directory instead of the ";
	print "default of oldTests/\n";
	print "\t--memory mem\n";
	print "\t\tIgnore files larger than mem\n";
	print "\t-f\n";
	print "\t\tPrint info only about failed tests\n";
	print Ci::helpFor("-nosu2");
	print Ci::helpFor("-h");
	exit(0);
}

defined($memory) or $memory = 50000000;
defined($failed) or $failed = 0;
defined($noSu2) or $noSu2 = 0;
defined($workdir) or $workdir = "tests";
defined($golddir) or $golddir = "oldTests";

my @tests = Ci::getTests();
my $total = scalar(@tests);

my @inRange = Ci::procRanges($ranges, $total);
my $rangesTotal = scalar(@inRange);

die "$0: No tests specified under -n\n" if ($rangesTotal == 0);

for (my $j = 0; $j < $rangesTotal; ++$j) {
        my $i = $inRange[$j];
        die "$0: out of range $i >= $total\n" if ($i >= $total);
	my $n = $tests[$i];

	my $isSu2 = Ci::isSu2("inputs/input$n.inp",$n);
	if ($isSu2 and $noSu2) {
		print STDERR "$0: WARNING: Ignored test $n ";
		print STDERR "because it's NOT an SU(2) test and ";
		print STDERR "you specified -nosu2\n";
		next;
        }

	procTest($n,$workdir,$golddir);
	print "-----------------------------------------------\n";
}

sub procTest
{
	my ($n,$workdir,$golddir) = @_;
	my %newValues;
	my %oldValues;
	procCout(\%newValues, $n,$workdir);
	procCout(\%oldValues, $n, $golddir);
	compareValues(\%newValues, \%oldValues, $n);
	procMemcheck($n);
}

sub procCout
{
	my ($values, $n, $dir) = @_;
	my $file = "$dir/runForinput$n.cout";
	if (!(-r "$file")) {
		print "|$n|: No $file found\n";
		return;
	}

	open(FILE, "$file") or return;
	my @energies;
	my $counter = 0;
	while (<FILE>) {
		chomp;
		if (/DMRG\+\+ version (.*$)/) {
			$values->{"version"} = $1;
			next;
		}

		if (/lowest eigenvalue= ([^ ]+) /) {
			push(@energies, $1);
			next;
		}
	}

	close(FILE);
	$values->{"energies"} = \@energies;
	$values->{"version"} = "UNDEFINED" unless (defined($values->{"version"}));

}

sub compareValues
{
	my ($newValues, $oldValues, $n) = @_;
	my $v1 = $newValues->{"version"};
	my $v2 = $oldValues->{"version"};
	defined($v1) or $v1 = "UNDEFINED";
	defined($v2) or $v2 = "UNDEFINED";
	print "|$n|: New Version $v1, Old Version $v2\n";
	my $maxEdiff = maxEnergyDiff($newValues->{"energies"}, $oldValues->{"energies"});
	print "|$n|: MaxEnergyDiff = $maxEdiff\n";
}

sub maxEnergyDiff
{
	my ($eNew, $eOld) = @_;
	return "NEW ENERGIES UNDEFINED" if (!defined($eNew));
	return "OLD ENERGIES UNDEFINED" if (!defined($eOld));
	my $n = scalar(@$eNew);
	return "ENERGY SIZES DIFFERENT" if (scalar(@$eOld) != $n);
	return "NO ENERGIES!" if ($n == 0);
	my $maxEdiff = 0;
	for (my $i = 0; $i < $n; ++$i) {
		my $tmp = abs($eNew->[$i] - $eOld->[$i]);
		$maxEdiff = $tmp if ($tmp > $maxEdiff);
 	}

	return "$maxEdiff [out of $n]";
}


sub procMemcheck
{
	my ($n) = @_;
	my $file1 = "$workdir/output$n.txt";
	my $mode = "OK";
	my $extra = "UNDEFINED";
	my %lost;
	open(FILE,$file1) or die "$0: Cannot open $file1 : $!\n";
	while (<FILE>) {
		if (/FATAL: You are requesting/ || /mandatory label/) {
			$mode = "FATAL";
			$extra = $_;
			chomp($extra);
			last;
		}

		if (/terminate called after throwing/) {
			$mode = "throw";
			while (<FILE>) {
				if (/what\(\)/) {
					$extra = $_;
					chomp($extra);
					last;
				}
			}

			last;
		}
		next unless (/^==/);
		if (/invalid/i) {
			$mode = "invalid";
			last;
		}

		if (/uninitial/i) {
			$mode = "uninitialized";
			last;
		}
		
		if (/([^ ]+) lost: ([^ ]+) bytes/) {
			my $val = $2;
			my $name = $1;
			$val =~ s/,//;
			$lost{"$name"}= $val;
			next;
		}
	}

	close(FILE);

	if ($mode eq "FATAL" || $mode eq "throw") {
		print "$0: ATTENTION TEST $n couldn't run because of $extra\n";
		return;
	}
	
	foreach my $key (keys %lost) {
		my $val = $lost{"$key"};
		next if ($val == 0);
		print "$0: ATTENTION TEST $n $key lost ".$lost{"$key"}." bytes\n";
	}

	return if ($mode eq "OK" and $failed);

	print "|$n|: output mode $mode\n";
}

sub fileSize
{
	my ($filename) = @_;
	my ($dev,$ino,$mode,$nlink,$uid,$gid,$rdev,$size,
	          $atime,$mtime,$ctime,$blksize,$blocks)
	              = stat($filename);
	return $size;
}

