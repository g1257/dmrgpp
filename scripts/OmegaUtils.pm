#!/usr/bin/perl

use strict;
use warnings;
use Math::Trig;

package OmegaUtils;

my $pi = Math::Trig::pi;

sub getLabels
{
	my ($hptr,$file) = @_;

	open(FILE,$file) or die "$0: Cannot open $file : $!\n";
	while (<FILE>) {
		chomp;
		foreach my $key (keys %$hptr) {
			if (/$key[= ]([^ ]+)/) {
				${$hptr->{$key}} = $1;
			}
		}
	}

	close(FILE);

	foreach my $key (keys %$hptr) {
		my $x = ${$hptr->{$key}};
		defined($x) or die "$0: Could not find $key in $file\n";
	}
}

sub printGnuplot
{
	my ($ptr, $geometry, $isPeriodic, $zeroAtCenter) = @_;

	my $factor = 0;
	my @fileIndices=(0);
	if ($geometry eq "chain") {
		$factor = 0.5;
	} elsif ($geometry eq "ladder") {
		$factor = 0.25;
		@fileIndices=(0,1);
	} else {
		die "$0: Unknown geometry $geometry\n";
	}

	foreach my $fileIndex (@fileIndices) {
		my $outFile = "outSpectrum$fileIndex.gnuplot";
		open(FOUT, ">", "$outFile") or die "$0: Cannot write to $outFile : $!\n";

		for my $omega (sort {$a <=> $b} keys %$ptr) {
			my $aptr = $ptr->{$omega};
			my $nks = scalar(@$aptr) - 1;
			my $numberOfQs = int($factor*$nks);
			my $centerShift = ($numberOfQs & 1) ? ($numberOfQs - 1)/2 : $numberOfQs/2;
			$centerShift = 0 unless ($zeroAtCenter);
			for (my $m2 = 0; $m2 < $numberOfQs; ++$m2) {
				my $m = $m2 - $centerShift;
				$m += $numberOfQs if ($m < 0);

				my $q = getQ($m2 - $centerShift, $numberOfQs, $isPeriodic);
				my $realPart = $aptr->[2*$m+1+2*$fileIndex*$numberOfQs];
				my $imagPart = $aptr->[2*$m+2+2*$fileIndex*$numberOfQs];
				print FOUT "$q $omega $realPart $imagPart\n";
			}
		}

		close(FOUT);
		print "$0: Written $outFile\n";
	}
}

sub printOffsetPlots
{
	my ($ptr, $geometry, $isPeriodic, $zeroAtCenter) = @_;

	my $factor = 0;
	my @fileIndices=(0);
	if ($geometry eq "chain") {
		$factor = 0.5;
	} elsif ($geometry eq "ladder") {
		$factor = 0.25;
		@fileIndices=(0,1);
	} else {
		die "$0: Unknown geometry $geometry\n";
	}

	foreach my $fileIndex (@fileIndices) {
		my $offset = 0.7*findMaxVertical($ptr, $factor, $fileIndex);
		my $outFile = "outSpectrum$fileIndex.offset";
		open(FOUT, ">", "$outFile") or die "$0: Cannot write to $outFile : $!\n";

		for my $omega (sort {$a <=> $b} keys %$ptr) {
			my $aptr = $ptr->{$omega};
			my $nks = scalar(@$aptr) - 1;
			my $numberOfQs = int($factor*$nks);
			my $centerShift = ($numberOfQs & 1) ? ($numberOfQs - 1)/2 : $numberOfQs/2;
			$centerShift = 0 unless ($zeroAtCenter);
			print FOUT "$omega ";
			for (my $m2 = 0; $m2 < $numberOfQs; ++$m2) {
				my $m = $m2 - $centerShift;
				$m += $numberOfQs if ($m < 0);

				my $q = getQ($m2 - $centerShift, $numberOfQs, $isPeriodic);
				#my $realPart = $aptr->[2*$m+1+2*$fileIndex*$numberOfQs];
				my $imagPart = $aptr->[2*$m+2+2*$fileIndex*$numberOfQs];
				my $val = $imagPart + $m2*$offset;
				print FOUT "$val ";
			}

			print FOUT "\n";
		}

		close(FOUT);
		print "$0: Written $outFile\n";
	}
}

sub findMaxVertical
{
	my ($ptr, $factor, $fileIndex) = @_;

	my $max = 0;
	for my $omega (sort {$a <=> $b} keys %$ptr) {
		my $aptr = $ptr->{$omega};
		my $nks = scalar(@$aptr) - 1;
		my $numberOfQs = int($factor*$nks);
		for (my $m = 0; $m < $numberOfQs; ++$m) {
			my $imagPart = abs($aptr->[2*$m + 2 + 2*$fileIndex*$numberOfQs]);
			$max = $imagPart if ($max < $imagPart);
		}
	}

	return $max;
}

sub getQ
{
	my ($m, $n, $isPeriodic) = @_;
	return ($isPeriodic) ? 2.0*$pi*$m/$n : ($m + 1)*$pi/($n+1.0);
}

1;

