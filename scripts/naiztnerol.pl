#!/usr/bin/perl

use strict;
use warnings;
use utf8;

#The file to be read has format
#omega value
#
my ($file, $delta) = @ARGV;
defined($file) or die "USAGE: $0 filename [delta]\n";
defined($delta) or $delta = 0.1;

print "#$file \n";
my @data = loadData($file);

my @data2 = findMaxima(\@data);

printMaxima(\@data2);

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
	my ($maxVal, $maxOmega, $prevVal) = (0, 0, 0);
	for (my $i = 0; $i < $n; ++$i) {
		my $ptr = $data->[$i];
		my $omega = $ptr->[0];
		my $val = $ptr->[1];

		if ($val > $prevVal) {
			#going up
			if ($val > $maxVal) {
				$maxVal = $val;
				$maxOmega = $omega;
			}
		} else {
			#going down
			if ($maxVal > 0) {
				push @data2, [$maxOmega, $maxVal];
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


