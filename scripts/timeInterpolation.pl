#!/usr/bin/env perl
#===============================================================================
#
#         FILE: timeEvolution3.pl
#
#        USAGE: ./timeEvolution3.pl
#
#  DESCRIPTION:
#
#      OPTIONS: ---
# REQUIREMENTS: ---
#         BUGS: ---
#        NOTES: ---
#       AUTHOR: G. A.
# ORGANIZATION:
#      VERSION: 1.0
#      CREATED: 08/14/2015
#     REVISION: ---
#===============================================================================

use strict;
use warnings;
use utf8;
my @m;
my @times;
my $numberOfSites;
my $timeIndex = 0;

while (<STDIN>) {
	chomp;
	my @temp = split;
	next if (/^#/);
	my $n = scalar(@temp);
	my $time = $temp[0];
	$times[$timeIndex] = $time;
	my $total = $n - 1;
	$total /= 2.0;
	my $index = 1;
	if (defined($numberOfSites)) {
		($total == $numberOfSites) or die "$0: Line $_\n";
	} else {
		$numberOfSites = $total;
	}

	for (my $i = 0; $i < $total; ++$i) {
		my @val = ($temp[$index],$temp[$index+1]);
		$index += 2;
		$m[$i + $timeIndex*$numberOfSites] = \@val;
	}

	$timeIndex++;
}

for (my $ti = 0; $ti < $timeIndex; ++$ti) {
	print "$times[$ti] ";
	for (my $i = 0; $i < $numberOfSites; ++$i) {
		my @val = getValue($i,$ti);
		my $n = scalar(@val);
		($n == 2) or die "$0: Problem with @val\n";
		my $b = ($val[0] eq "0.00" or $val[1] eq "0.00");
		if ($b and $ti + 1 < $timeIndex) {
			($ti > 0) or die "$0: Cannot interpolate for time=0\n";
			($ti + 1 < $timeIndex) or die "$0: Cannot interpolate for last time\n";
			my @val1 = getValue($i,$ti-1);
			my @val2 = getValue($i,$ti+1);
			@val = averageOf(\@val1,\@val2);
		}

		print "@val ";
	}

	print "\n";
}

sub getValue
{
	my ($i,$ti) = @_;
	my $ptr = $m[$i + $ti*$numberOfSites];
	return @$ptr;
}

sub averageOf
{
	my ($ptr1,$ptr2) = @_;
	my @val1 = @$ptr1;
	my @val2 = @$ptr2;
	my $n1 = scalar(@val1);
	my $n2 = scalar(@val2);
	($n1 == $n2) or die "$0: Cannot average arrays of different size\n";
	my @result;
	for (my $i = 0; $i < $n1; ++$i) {
		@result[$i] = 0.5*($val1[$i] + $val2[$i]);
	}

	return @result;
}

