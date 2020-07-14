#!/usr/bin/perl

use strict;
use warnings;
use utf8;

my ($file) = @ARGV;
defined($file) or die "USAGE: $0 filename\n";

my %data = loadData($file);
die "$0: No k values found in $file\n" if (scalar(keys %data) == 0);

printData(\%data);

sub printData
{
	my ($data) = @_;
	foreach my $k (sort {$a <=> $b} keys %$data) {
		print "$k ".$data->{$k}."\n";
	}
}

sub loadData
{
	my ($file) = @_;
	my %data;
	open(FILE, "<", $file) or die "$0: Cannot open $file : $!\n";
	while (<FILE>) {
		chomp;
		my @temp = split;
		my $n = scalar(@temp);
		next if ($n != 3);
		my $k = $temp[0];
		my $omega = $temp[1];
		my $value = $temp[2];
		my $sum = $data{$k};
		if (defined($sum)) {
			$sum += $value;
		} else {
			$sum = $value;
		}

		$data{$k} = $sum;
	}

	close(FILE);
	return %data;
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



