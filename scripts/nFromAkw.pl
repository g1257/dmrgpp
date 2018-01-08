#!/usr/bin/perl

use strict;
use warnings;
use utf8;
use Math::Trig;
use OmegaUtils;

my $pi = Math::Trig::pi;

use Getopt::Long qw(:config no_ignore_case);

my $usage = "-f dollarizedInput [-p] [-z] minusSpectrum plusSpectrum\n";

my $templateInput;
my $isPeriodic;
my $zeroAtCenter = 0;

GetOptions('f=s' => \$templateInput,
           'p' => \$isPeriodic,
           'z' => \$zeroAtCenter) or die "$usage\n";

(defined($templateInput) and defined($isPeriodic)) or die "$0: USAGE: $usage\n";

my $geometry;
my $hptr = {"GeometryKind" => \$geometry};
OmegaUtils::getLabels($hptr,$templateInput);

my ($fileMinus, $filePlus) = @ARGV;

my $numberKs;
my %specMinus;
readSpectrum(\%specMinus, \$numberKs, $fileMinus);

my %specPlus;
readSpectrum(\%specPlus, \$numberKs, $filePlus);

if ($geometry eq "ladder") {
	my @nkxpi;
	sumOverOmega(\@nkxpi, \%specMinus, 1);
	printVsQ("outnkxpi.dat", \@nkxpi);
}

my @nkx0;
sumOverOmega(\@nkx0, \%specMinus, 0);
printVsQ("outnkx0.dat", \@nkx0);

my $totalMy = ($geometry eq "ladder") ? 2 : 1;

for (my $mp = 0; $mp < 2; ++$mp) { #mp = 0 is -, mp=1 is +
	my $ptr = ($mp == 0) ? \%specMinus : \%specPlus;
	for (my $my = 0; $my < $totalMy; ++$my) {
		my %h;
		sumOverKx(\%h, $ptr, $my);
		printVsOmega("nVsOmegaky$my"."Sector$mp.dat", \%h);
	}
}

sub printVsOmega
{
	my ($fout, $ptr) = @_;
	open(FOUT, ">", "$fout") or die "$0: Cannot write to $fout : $!\n";
	for my $omega (sort {$a <=> $b} keys %$ptr) {
		my $val = $ptr->{$omega};
		$val = 0 if ($val < 0);
		print FOUT "$omega $val\n";
	}

	close(FOUT);
	print STDERR "$0: File $fout has been written.\n";
}

sub sumOverKx
{
	my ($v, $ptr, $my) = @_;
	my $factor = 0;
	my @fileIndices=(0);
	if ($geometry eq "chain") {
		$factor = 0.5;
		die "$0: Chain does not have ky != 0\n" if ($my != 0)
	} elsif ($geometry eq "ladder") {
		$factor = 0.25;
		@fileIndices=(0,1);
	} else {
		die "$0: Unknown geometry $geometry\n";
	}

	my $fileIndex = $my;

	for my $omega (sort {$a <=> $b} keys %$ptr) {
		my $aptr = $ptr->{$omega};
		my $nks = scalar(@$aptr) - 1;
		my $numberOfQs = int($factor*$nks);
		for (my $m2 = 0; $m2 < $numberOfQs; ++$m2) {
			my $realPart = $aptr->[2*$m2+1+2*$fileIndex*$numberOfQs];
			#my $imagPart = $aptr->[2*$m2+2+2*$fileIndex*$numberOfQs];
			if (defined($v->{$omega})) {
				$v->{$omega} += $realPart;
			} else {
				$v->{$omega} = $realPart;
			}
		}
	}
	
}

sub printVsQ
{
	my ($fout, $v) = @_;
	my $numberOfQs = scalar(@$v);
	my $centerShift = ($numberOfQs & 1) ? ($numberOfQs - 1)/2 : $numberOfQs/2;
	$centerShift = 0 unless ($zeroAtCenter);

	open(FOUT, ">", "$fout") or die "$0: Cannot write to $fout : $!\n";
	for (my $m2 = 0; $m2 < $numberOfQs; ++$m2) {
		my $m = $m2 - $centerShift;
		$m += $numberOfQs if ($m < 0);
		my $q = getQ($m2 - $centerShift, $numberOfQs, $isPeriodic);
		print FOUT "$q ".$v->[$m2]."\n";
	}

	close(FOUT);
	print STDERR "$0: File $fout has been written.\n";
}

sub sumOverOmega
{

	my ($v, $ptr, $my) = @_;
	
	my $factor = 0;
	my @fileIndices=(0);
	if ($geometry eq "chain") {
		$factor = 0.5;
		die "$0: Chain does not have ky != 0\n" if ($my != 0)
	} elsif ($geometry eq "ladder") {
		$factor = 0.25;
		@fileIndices=(0,1);
	} else {
		die "$0: Unknown geometry $geometry\n";
	}

	my $fileIndex = $my;

	for my $omega (sort {$a <=> $b} keys %$ptr) { #no need to sort
		my $aptr = $ptr->{$omega};
		my $nks = scalar(@$aptr) - 1;
		my $numberOfQs = int($factor*$nks);
		for (my $m2 = 0; $m2 < $numberOfQs; ++$m2) {
			my $realPart = $aptr->[2*$m2+1+2*$fileIndex*$numberOfQs];
			#my $imagPart = $aptr->[2*$m2+2+2*$fileIndex*$numberOfQs];
			if (defined($v->[$m2])) {
				$v->[$m2] += $realPart;
			} else {
				$v->[$m2] = $realPart;
			}
		}
	}
}

sub readSpectrum
{
	my ($ptr, $ptrN, $file) = @_;

	my $isGood = 1;
	open(FILE, "<", $file) or die "$0: Cannot open $file : $!\n";
	while (<FILE>) {
		next if (/^#/);
		my @temp = split;
		my $n = scalar(@temp);
		if ($n < 2) {
			print STDERR "$0: Ignored line $. in $file, less than 2 cols\n";
			next;
		}

		$$ptrN = $n - 1 if (!defined($$ptrN));
		if ($n - 1 != $$ptrN) {
			$isGood = 0;
			print STDERR "$0: Line $. in $file with wrong cols\n";
			last;
		}

		my $omega = $temp[0];
		my $oldVal = $ptr->{$omega};
		if (defined($oldVal)) {
			for (my $i = 1; $i < $$ptrN; ++$i) {
				$temp[$i] += $oldVal->[$i];
			}
		}

		$ptr->{$omega} = \@temp;
	}

	close(FILE);

	return if ($isGood);

	die "$0: $file with at least 1 line with wrong number of cols, ".$$ptrN." expected\n";
}

sub getQ
{
	my ($m, $n, $isPeriodic) = @_;
	return ($isPeriodic) ? 2.0*$pi*$m/$n : $m*$pi/($n+1.0);
}


