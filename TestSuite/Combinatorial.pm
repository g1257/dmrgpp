#!/usr/bin/perl

use strict;
use warnings;
use utf8;

package Combinatorial;

sub main
{
	my ($levels, $callback) = @_;
	my $total = 1;
	my @sizes;
	my $counter = scalar(@$levels);
	for (my $i = 0; $i < $counter; ++$i) {
		my $ptr = $levels->[$i];
		my $n = scalar(@$ptr);
		$total *= $n;
		push @sizes, $n;
	}

	print STDERR "$0: Total configurations $total\n";

	for (my $i = 0; $i < $total; ++$i) {
		my @location = findLocation($i, \@sizes, $total);
		my @items = findItemsForLocation(\@location, $levels);
		$callback->(\@items);
	}
}

sub findLocation
{
	my ($ind, $sizes, $volume) = @_;
	# ind = x[0] + x[1]*sizes[0] + x[2]*sizes[0]*sizes[1] + ...
	my @location;
	my $n = scalar(@$sizes);
	die "$0: findLocation sizes==0\n" if ($n == 0);

	my $factor = $volume;

	my $tmp = $ind;
	my $j = $n - 1;
	for (my $i = 0; $i < $n; ++$i) {
		$factor /= $sizes->[$j];
		$location[$j] = int($tmp/$factor);
		$tmp -= $location[$j--]*$factor;
	}

	return @location;
}

sub findItemsForLocation
{
	my ($location, $levels) = @_;

	my $n = scalar(@$location);
	if ($n != scalar(@$levels)) {
		die "$0: findItemsForLocation\n";
	}

	my @items;
	for (my $i = 0; $i < $n; ++$i) {
		my $j = $location->[$i];
		my $item = $levels->[$i]->[$j];
		push @items, $item;
	}

	return @items;
}

1;

