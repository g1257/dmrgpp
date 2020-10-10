#!/usr/bin/perl

use strict;
use warnings;
use utf8;

my @sum;
foreach my $file (@ARGV) {
	my @array;
	loadFile(\@array, $file);
	addFile(\@sum, \@array);
}

printTotal(\@sum);

sub printTotal
{
	my ($array) = @_;
	my $n = scalar(@$array);
	my $prevOmega = 0;
	for (my $i = 0; $i < $n; ++$i) {
		my $ptr = $array->[$i];
		my $omega = $ptr->[1];
		if ($i > 0 && abs($omega - $prevOmega)>1e-6) {
			print "\n";
		}

		print "@$ptr\n";
		$prevOmega = $omega;
	}
}

sub loadFile
{
	my ($h, $file) = @_;
	open(FILE, "<", $file) or die "$0: Cannot open file $file : $!\n";
	while (<FILE>) {
		chomp;
		next if ($_ eq "");
		my @temp = split;
		die "$0: Incorrect line $_ in $file\n" if (scalar(@temp) != 3);
		push @$h, \@temp;
	}

	close(FILE);
}

sub addFile
{
	my ($sum, $array) = @_;
	my $n = scalar(@$array);
	my $m = scalar(@$sum);
	if ($m == 0) {
		@$sum = @$array;
		return;
	}

	if ($m != $n) {
		die "$0: Trying to add an array of different size to previous $n != $m\n";
	}

	for (my $i = 0; $i < $n; ++$i) {
		my $ptr1 = $array->[$i];
		my $ptr2 = $sum->[$i];
		my $nn = scalar(@$ptr1);
		die "$0: addFile: $i-th array incorrect size\n" if ($nn != scalar(@$ptr2));
		die "$0: addFile: Neither $i-th array is size 3\n" if ($nn != 3);
		my $k = $ptr1->[0];
		my $diff0 = abs($k - $ptr2->[0]);
		die "$0: addFile: Wrong k point\n" if ($diff0 > 1e-5);
		my $omega = $ptr1->[1];
		my $diff1 = abs($omega - $ptr2->[1]);
		die "$0: addFile: Wrong omega point\n" if ($diff1 > 1e-5);
		$ptr2->[2] += $ptr1->[2];
	}
}

