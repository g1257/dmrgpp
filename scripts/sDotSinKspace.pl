#!/usr/bin/perl

use strict;
use warnings;
use utf8;
use Fourier;

my ($filename) = @ARGV;
defined($filename) or die "USAGE: $0 filename\n";

my $labelSzSz = "<gs|sz;sz|gs>";
my @szsz = Fourier::loadMatrix($filename, $labelSzSz);

my $labelSpSm = "<gs|splus;sminus|gs>";
my @spsm = Fourier::loadMatrix($filename, $labelSpSm);

my @sDotS = createSdotS(\@szsz, \@spsm);

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
	my ($szsz, $spsm) = @_;
	my @sDotS;
	my $rows = scalar(@$szsz);
	if ($rows != scalar(@spsm)) {
		die "$0: Sz.Sz and Sp.Sm arrays have different row sizes\n";
	}

	for (my $i = 0; $i < $rows; ++$i) {
		my $szszPointer = $szsz->[$i];
		my $spsmPointer = $spsm->[$i];
		my $cols = scalar(@$szszPointer);
		if ($cols != scalar(@$spsmPointer)) {
			die "$0: Sz.Sz and Sp.Sm arrays have different row cols for row $i\n";
		}

		my @temp;
		for (my $j = 0; $j < $cols; ++$j) {
			$temp[$j] = $szszPointer->[$j] + $spsmPointer->[$j];
		}

		$sDotS[$i] = \@temp;
	}

	return @sDotS;
}

