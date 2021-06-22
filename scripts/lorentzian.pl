#!/usr/bin/perl

use strict;
use warnings;
use utf8;

#The file to be read has format
#weight0 omega0
#weight1 omega1
#...
#
my ($file, $wbegin, $wtotal, $wstep, $delta, $norm) = @ARGV;
defined($delta) or die "USAGE: $0 filename omegaBegin omegaTotal omegaStep delta [norm]\n";

print "#$file, $wbegin, $wtotal, $wstep, $delta\n";
my @data = loadData($file);

my $wstruct = {
	"begin" => $wbegin,
	"total" => $wtotal,
	"step" => $wstep,
	"delta" => $delta};

my @data2 = computeData(\@data, $wstruct);

my $sumptr = pop @data2;

my $sum = $sumptr->[1];

my $factor = defined($norm) ? $norm/$sum : 1;

plotData(\@data2, $factor);

sub plotData
{

	my ($data, $factor) = @_;
	my $n = scalar(@$data);
	for (my $i = 0; $i < $n; ++$i) {
		my $ptr = $data->[$i];
		my $omega = $ptr->[0];
		my $val = $factor * $ptr->[1];
		print "$omega $val\n";
	}
}

sub computeData
{
	my ($data, $wstruct) = @_;
	my @data2;
	my $sum = 0;
	for (my $i = 0; $i < $wstruct->{"total"}; ++$i) {
		my $omega = $wstruct->{"begin"} + $wstruct->{"step"}*$i;
		my $val = computeValueFor($omega, $data, $wstruct->{"delta"});
		$data2[$i] = [$omega, $val];
		$sum += $val;
		#print "$omega $val\n";
	}

	push @data2, [0, $sum];
	return @data2;
}

sub computeValueFor
{
	my ($omega, $data, $delta) = @_;

	my $sum = 0;
	my $total = scalar(@$data);
	for (my $i = 0; $i < $total; ++$i) {
		my $ptr = $data->[$i];
		my $deltaOmega = $omega - $ptr->[1];
		$sum += $ptr->[0]/($deltaOmega*$deltaOmega + $delta*$delta);
	}

	return $sum;
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

	print STDERR "$0: Found $count peaks\n";
	return @data;
}


