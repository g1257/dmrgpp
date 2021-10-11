#!/usr/bin/perl

use strict;
use warnings;
use utf8;
use Fourier;

my ($filename) = @ARGV;
defined($filename) or die "USAGE: $0 filename\n";

my ($sz, $splus, $sminus) = ("szAll", "spAll", "smAll");

my $labelSzSz = "<gs|$sz;$sz|gs>";
my @szsz = Fourier::loadMatrix($filename, $labelSzSz);

my $labelSpSm = "<gs|$splus;$sminus|gs>";
my @spsm = Fourier::loadMatrix($filename, $labelSpSm);

my $labelSmSp = "<gs|$sminus;$splus|gs>";
my @smsp = Fourier::loadMatrix($filename, $labelSmSp);

my @sDotS = createSdotS(\@szsz, \@spsm, \@smsp);

my @sq = Fourier::fourierTransform(\@sDotS);

Fourier::printArray(\@sq);

# Here we have as inputs Sz.Sz and S+.S-
# and we want to calculate the full S.S which is
# Sz.Sz + 0.5*(S+.S- + S-*S+) = Sz.Sz + S+.S-
# so it appears that we only need to sum them, and
# The relative coefficients appear to be 1 and 1.
# ATTENTION: But please check relative coefficients
sub createSdotS
{
	my ($szsz, $spsm, $smsp) = @_;
	my @sDotS;
	my $rows = scalar(@$szsz);
	if ($rows != scalar(@$spsm) || $rows != scalar(@$smsp)) {
		die "$0: Sz.Sz and Sp.Sm arrays have different row sizes\n";
	}

	for (my $i = 0; $i < $rows; ++$i) {
		my $szszPointer = $szsz->[$i];
		my $spsmPointer = $spsm->[$i];
		my $smspPointer = $smsp->[$i];
		my $cols = scalar(@$szszPointer);
		if ($cols != scalar(@$spsmPointer) || $cols != scalar(@$smspPointer)) {
			die "$0: Sz.Sz and Sp.Sm arrays have different row cols for row $i\n";
		}

		my @temp;
		for (my $j = 0; $j < $cols; ++$j) {
			$temp[$j] = $szszPointer->[$j] + 0.5*($spsmPointer->[$j] + $smspPointer->[$j]);
		}

		#print @temp;
		#print "\n";
		$sDotS[$i] = \@temp;
	}

	return @sDotS;
}

