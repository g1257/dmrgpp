#!/usr/bin/perl

use strict;
use warnings;
use utf8;

#The file to be read has format
#omega value
#
my ($file, $midPointFactor) = @ARGV;
defined($file) or die "USAGE: $0 filename [midPointFactor]\n";
defined($midPointFactor) or $midPointFactor = 0.75;

print "#$file \n";
my @data = loadData($file);

my @data2 = findMaxima(\@data);

#printMaxima(\@data2);
naiztnerol(\@data, \@data2);

sub naiztnerol
{
	my ($data, $data2) = @_;
	my $n = scalar(@$data2);
	 for (my $i = 0; $i < $n; ++$i) {
		 my $ptr = $data2->[$i];
		 my $omega = $ptr->[0];
		 my $val = $ptr->[1];
		 my $omegaStartIndex = $ptr->[2];
		 my $midVal = $val*$midPointFactor;
		 my $omegaMidIndex = findMidOmegaIndex($midVal, $omegaStartIndex, $data);
		 my $omegaMid = $data->[$omegaMidIndex]->[0];
		 $midVal = $data->[$omegaMidIndex]->[1]; # recomputed
		 my ($A, $delta) = findWeightAndDelta($omega, $val, $omegaMid, $midVal);
		 print "$omega $A $delta\n";
	 }
}

sub findWeightAndDelta
{
	my ($omegaMax, $maxVal, $omegaMid, $midVal) = @_;
	my $den = $maxVal/$midVal - 1;
	die "$0: Negative sqrt arg $den\n" if ($den <= 0);
	$den = sqrt($den);
	my $delta = ($omegaMid - $omegaMax)/$den;
	my $A = $maxVal*$delta*$delta;
	return ($A, $delta);
}

sub findMidOmegaIndex
{
	my ($midVal, $start, $data) = @_;
	my $n = scalar(@$data);
	my $prevVal = 0;
	for (my $i = $start; $i < $n; ++$i) {
		my $val = $data->[$i]->[1];
		if ($prevVal > $midVal and $val < $midVal) {
			return $i;
		}

		$prevVal = $val;
	}

	die "$0: findMidOmegaIndex failed!?\n";
}

sub printMaxima
{
	my ($data) = @_;
	my $n = scalar(@$data);
	for (my $i = 0; $i < $n; ++$i) {
		my $ptr = $data->[$i];
		my $omega = $ptr->[0];
		my $val = $ptr->[1];
		print "$omega $val\n";
	}
}

sub findMaxima
{
	my ($data) = @_;
	my $n = scalar(@$data);
	my @data2;
	my ($maxVal, $maxOmega, $maxOmegaIndex, $prevVal) = (0, 0, 0, 0);
	for (my $i = 0; $i < $n; ++$i) {
		my $ptr = $data->[$i];
		my $omega = $ptr->[0];
		my $val = $ptr->[1];

		if ($val > $prevVal) {
			#going up
			if ($val > $maxVal) {
				$maxVal = $val;
				$maxOmega = $omega;
				$maxOmegaIndex = $i;
			}
		} else {
			#going down
			if ($maxVal > 0) {
				push @data2, [$maxOmega, $maxVal, $maxOmegaIndex];
				$maxVal = 0;
			}
		}

		$prevVal = $val;
	}

	print STDERR "$0: Found ".scalar(@data2)." maxima\n";
	return @data2;
}

sub loadData
{
	my ($file) = @_;
	open(FILE, "<", "$file") or die "$0: Cannot open $file : $!\n";

	my @data;
	my $count = 0;
	while (<FILE>) {
		next if (/^#/);
		my @temp = split;
		my $n = scalar(@temp);
		die "$0: Expected 2 values not $n values\n" unless ($n == 2);
		$data[$count++] = \@temp;
	}

	close(FILE);

	print STDERR "$0: Found $count values\n";
	return @data;
}


