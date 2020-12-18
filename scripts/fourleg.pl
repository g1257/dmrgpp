#!/usr/bin/perl

use strict;
use warnings;
use utf8;

use OmegaUtils;

my ($file, $l) = @ARGV;
defined($l) or die "USAGE: $0 filename l\n";

my @array = loadValues($file);
ft(\@array, $l);

sub ft
{
	my ($a, $l) = @_;
	my ($isPeriodic, $zeroAtCenter, $nonNegativeOnly) = (1, 0, 0);
	my $omegas = scalar(@$a);
	my $centralSite = $l;
	my $geometry = {"name" => "ladder", "leg" => 2, "subname" => ""};
	my $hptr = {isPeriodic => $isPeriodic, centralSite => $centralSite};
	my $outSpectrum = "out.spectrum";

	open(FOUTSPECTRUM, ">", "$outSpectrum") or die "$0: Cannot write to $outSpectrum : $!\n";

	for (my $index = 0; $index < $omegas; ++$index) {
		my $ptr = $a->[$index];
		(scalar(@$ptr) == 2) or die "$0: ptr.size should be 2\n";
		my $omega = $ptr->[0];
		print FOUTSPECTRUM "$omega ";
		my $ptr2 = $ptr->[1];
		my @v;
		my $sites = scalar(@$ptr2);
		my $counter = 0;
		for (my $i = 0; $i < $sites; ++$i) {
			my $y = $i % 4;
			next if ($y > 1);
			$v[$counter++] = $ptr2->[$i];
		}

		die "$0: Counter error\n" if ($counter != 2*$l);

		my @f;
		OmegaUtils::fourierLadder(\@f, \@v, 2, $hptr);
		my @array;
		OmegaUtils::writeFourier(\@array, \@f, $geometry);
		printSpectrum(\@array);
	}

	close(FOUTSPECTRUM);
	OmegaUtils::printGnuplot($outSpectrum, $geometry, $isPeriodic, $zeroAtCenter, $nonNegativeOnly);
}

sub printSpectrum
{
	my ($array) = @_;

	my $n = scalar(@$array);
	for (my $j = 0; $j < $n; ++$j) {
		my $array2 = $array->[$j];
		my @array2 = @$array2;
		print FOUTSPECTRUM "$array2[1] $array2[2] ";
	}

	print FOUTSPECTRUM "\n";
}

sub loadValues
{
	my ($file) = @_;
	open(FILE, "<", $file) or die "$0: Cannot open $file : $!\n";
	my @array;
	my $index = 0;
	while (<FILE>) {
		my @temp = split;
		my $n = scalar(@temp);
		die "$0: Expected two values, got $n\n" unless ($n == 2);
		my $omega = $temp[0];
		my $sites = $temp[1];
		my @values;
		for (my $i = 0; $i < $sites; ++$i) {
			$_ = <FILE>;
			chomp;
			my @temp = split;
			my $x = scalar(@temp);
			die "$0: Expected three values, got $x\n" unless ($x == 3);
			die "$0: Expected $i as first number, got $temp[0]\n" unless ($i == $temp[0]);
			my @temp2 = ($temp[1], $temp[2]);
			$values[$i] = \@temp2;
		}

		my @tmp = ($omega, \@values);
		$array[$index++] = \@tmp;
	}

	close(FILE);
	print STDERR "$0: Found $index omegas\n";
	return @array;
}

