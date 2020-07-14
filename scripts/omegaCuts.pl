#!/usr/bin/perl

use strict;
use warnings;
use utf8;

my ($file) = @ARGV;
defined($file) or die "USAGE: $0 filename\n";

my @kvalues = (0, 3.1415927);
my @data = loadData($file, \@kvalues);
die "$0: No k values found in $file\n" if (scalar(@data) == 0);

printData(\@data, \@kvalues);

sub printData
{
	my ($data) = @_;
	my $n = scalar(@$data);
	die "$0: No k values found\n" if ($n == 0);
	print STDERR "$0: Found $n kvalues\n";
	my $h = $data->[0];
	foreach my $omega (sort {$a <=> $b} keys %$h) {
		print "$omega";
		for (my $j = 0; $j < $n; ++$j) {
			my $value = $data->[$j]->{$omega};
			die "$0: Value too negative\n" if ($value < -4);
			$value = 0 if ($value < 0);
			print " $value";
		}

		print "\n";
	}
}

sub loadData
{
	my ($file, $kvalues) = @_;
	my @data;
	open(FILE, "<", $file) or die "$0: Cannot open $file : $!\n";
	while (<FILE>) {
		chomp;
		my @temp = split;
		my $n = scalar(@temp);
		next if ($n != 3);
		my $k = $temp[0];
		my $omega = $temp[1];
		my $value = $temp[2];
		my $ind = findKindex($kvalues, $k);
		next if ($ind < 0);
		my $hptr = $data[$ind];
		if (defined($hptr)) {
			my $sum = $hptr->{$omega};
			if (defined($sum)) {
				$sum += $value;
			} else {
				$sum = $value;
			}

			$hptr->{$omega} = $sum;
		} else {
			my $h = {"$omega" => $value };
			$data[$ind] = $h;
		}
	}

	close(FILE);
	return @data;
}

sub findKindex
{
	my ($kvalues, $k) = @_;
	my $n = scalar(@$kvalues);
	for (my $i = 0; $i < $n; ++$i) {
		return $i if (abs($k - $kvalues->[$i]) < 1e-5);
	}
	
	return -1;
}



