#!/usr/bin/perl

use strict;
use warnings;

my ($omega) = @ARGV;
defined($omega) or die "USAGE: $0 omega\n";

my $maxSite = 0;
my @values;
my @values2;
my $status;

while(<STDIN>) {
	if (/P3/) {
		$status="p3";
	} elsif (/P2/) {
		$status="p2";
	} else {
		next;
	}

	chomp;
	my @temp = split;
	die "$0: Line $_\n" unless (scalar(@temp)==5);

	my $site = $temp[0];
	$values[$site] = $temp[1] if ($status eq "p3");
	$values2[$site] = $temp[1] if ($status eq "p2");
	$maxSite = $site if ($maxSite < $site);
}

$maxSite++;
print STDERR "$0: Read $maxSite sites\n";

print "#omega=$omega\n";
for (my $i = 0; $i < $maxSite; ++$i) {
	my $v = $values[$i];
	my $v2 = $values2[$i];
	defined($v) or die "$0: Undefined value for site = $i\n";
	print "$i $v $v2\n";
}

