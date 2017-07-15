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

sub helpFor
{
	my ($label) = @_;
	my $h = "";
	if ($label eq "-n") {
		$h .= "\t-n n\n";
		$h .= "\t\tIgnore all tests except test(s) supplied\n";
		$h .= "\t\tThis is a comma-separated list of at least one range.\n";
		$h .= "\t\tA range is one of the following.\n";
		$h .= "\t\t\tA number, like 2\n";
		$h .= "\t\t\tA number followed by a dash, like 2-; this sets the minimum\n";
		$h .= "\t\t\tA dash followed by a number, like -2; this sets the maximum\n";
		$h .= "\t\t\tTwo numbers separated by a dash, like 2-4, indicating the range {2, 3, 4}\n";
		return $h;
	} elsif ($label eq "-w") {
		$h .= "\t-w workdir\n";
		$h .= "\t\tUse workdir as working directory not the default of tests/\n";
		return $h;
	} elsif ($label eq "-nosu2") {
		$h .= "\t-nosu2\n";
		$h .= "\t\tDo not postprocess SU(2) tests\n";
		return $h;
	} elsif ($label eq "-h") {
		$h .= "\t-h\n";
		$h .= "\t\tPrint this help and exit\n";
		return $h;
	}

	die "$0: No printHelpFor $label\n";
}

1;

