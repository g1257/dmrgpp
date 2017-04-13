#!/usr/bin/perl

use strict;
use warnings;
use Getopt::Long qw(:config no_ignore_case);
use Ci;

my ($min,$max,$memory,$failed,$help);
my ($workdir,$golddir,$n);
GetOptions(
'm=i' => \$min,
'M=i' => \$max,
'n=i' => \$n,
'memory=i' => \$memory,
'f' => \$failed,
'h' => \$help,
'g=s' => \$golddir,
'w=s' => \$workdir);

if (defined($help)) {
	print "USAGE: $0 [options]\n";
	print "\tIf no option is given examines all runs\n";
	print "\t-m min\n";
	print "\t\tMinimum test to examine is min (inclusive)\n";
	print "\t-M max\n";
	print "\t\tMaximum test to examine is max (inclusive)\n";
	print "\t-n n\n";
	print "\t\tIgnore all tests except test number n\n";
	print "\t-w workdir\n";
	print "\t\tUse workdir as working directory not the default of tests/\n";
	print "\t-g golddir\n";
	print "\t\tUse golddir for oracles directory instead of the ";
	print "default of oldTests/\n";
	print "\t--memory mem\n";
	print "\t\tIgnore files larger than mem\n";
	print "\t-f\n";
	print "\t\tPrint info only about failed tests\n";
	print "\t-h\n";
	print "\t\tPrint this help and exit\n";
	exit(0);
}

defined($memory) or $memory = 50000000;
defined($failed) or $failed = 0;
defined($workdir) or $workdir = "tests";
defined($golddir) or $golddir = "oldTests";
if (defined($n)) {
	if (defined($min) or defined($max)) {
		die "$0: -n cannot be used with either -m or -M\n";
	}

	$min = $n;
	$max = $n;
}

my @tests;
Ci::getTests(\@tests);

my $total = scalar(@tests);
for (my $i = 0; $i < $total; ++$i) {
	my $n = $tests[$i];
	next if (defined($min) and $n < $min);
	next if (defined($max) and $n > $max);
	procTest($n,$workdir,$golddir);
	print "-----------------------------------------------\n";
}

sub procTest
{
	my ($n,$workdir,$golddir) = @_;
	procData($n,$workdir,$golddir);
	procMemcheck($n);
}

sub procData
{
	my ($n,$workdir,$golddir) = @_;
	my $file1 = "$workdir/data$n.txt";
	my $file2 = "$golddir/data$n.txt";
	(-r "$file1") or return;
	(-r "$file2") or return;
	my $size1 = fileSize($file1);
	my $size2 = fileSize($file2);
	$size1 = $size2 if ($size2 > $size1);
	if ($memory > 0 and $size2 > $memory) {
		print STDERR "$0: File $file1 or $file2 is too big (ignoring test)\n";
		return;
	}

	my $cmd = "diff $file1 $file2";
	my @version;
	my ($newEnergy, $oldEnergy);
	my $maxEdiff = 0;
	open(PIPE,"$cmd |") or return;
	while (<PIPE>) {
		chomp;
		if (/([\<\>]) DMRG\+\+ version (.*$)/) {
			my $tmp = ($1 eq "<") ? 0 : 1;
			$version[$tmp] = $2;
			next if (!defined($version[0]) or !defined($version[1]));
			print "|$n|: New Version $version[0], Old Version $version[1]\n";
			next;
		}

		if (/\< \#Energy=(.+$)/) {
			$newEnergy = $1;
			next;
		}

		if (/\> \#Energy=(.+$)/ and defined($newEnergy)) {
			$oldEnergy = $1;
			my $tmp = $newEnergy-$oldEnergy;
			print "|$n|: EnergyNew-EnergyOld=$tmp\n";
			$tmp = abs($tmp);
			$maxEdiff = $tmp if ($maxEdiff < $tmp);
			undef($newEnergy);
			next;
		}

	}

	close(PIPE);
	print "|$n|: MaxEnergyDiff = $maxEdiff\n" if ($maxEdiff > 0);
}

sub procMemcheck
{
	my ($n) = @_;
	my $file1 = "tests/output$n.txt";
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

