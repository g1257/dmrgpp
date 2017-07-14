#!/usr/bin/perl

use strict;
use warnings;

package Ci;

sub getTests
{
	my $desc = "inputs/descriptions.txt";
	my @tests;
	my $counter = 0;
	open(FILE,"$desc") or die "$0: Cannot open $desc : $!\n";

	while (<FILE>) {
		last if (/^\#TAGSTART/);
	}

	while (<FILE>) {
		last if (/^\#TAGEND/);
		next if (/^\#/);
		if (/^(\d+)\)/) {
			$tests[$counter++] = $1;
		}
	}

	close(FILE);
	return @tests;
}

sub isSu2
{
	my ($file,$n) = @_; 
	open(FILE, "$file") or return 0;
	my $su2 = 0;
	while (<FILE>) {
		chomp;
		if (/UseSu2Symmetry=1/) {
			$su2 = 1;
			last;
		}
	}

	close(FILE);
	return $su2;
}

sub procRanges
{
	my ($range, $total) = @_;
	my @inRange;
	if (!defined($range)) {
		for (my $i = 0; $i < $total; ++$i) {
			push @inRange, $i;
		}

		return @inRange;
	}

	my @temp = split(/,/, $range);
	my $n = scalar(@temp);
	for (my $i = 0; $i < $n; ++$i) {
		procRange(\@inRange, $temp[$i], $total);
	}

	my @unique = do { my %seen; grep { !$seen{$_}++ } @inRange };
	return @unique;
}

sub procRange
{
	my ($ranges, $range, $total) = @_;
	my @temp = split(/\-/, $range);
	my $n = scalar(@temp);
	die "$0: FATAL: Empty range $range\n" if ($n == 0);
	if ($n == 1 and $range =~ /^[0-9]+$/) {
		push @$ranges, $range;
		return;
	}

	if ($n == 2 and $temp[0] =~ /^[0-9]+$/ and $temp[1] =~ /^[0-9]+$/) {
		for (my $i = $temp[0]; $i <= $temp[1]; ++$i) {
			push @$ranges, $i;
		}

		return;
	}

	die "$0: FATAL: Invalid range spec.: $range\n";
}

1;

