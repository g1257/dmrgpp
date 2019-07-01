#!/usr/bin/perl

use strict;
use warnings;
use utf8;
use Math::Trig;
use lib ".";
use OmegaUtils;

my $pi = Math::Trig::pi;

use Getopt::Long qw(:config no_ignore_case);

=pod
Place the out.spectrum for minus spectrum into (say) outMinus.spectrum
and   the out.spectrum for plus  spectrum into (say) outPlus.spectrum
Compute the mu = E(Nup, Ndown) - E(Nup - 1, Ndown) or with the +1 particle
Take one input dollarized file (either one is fine) and run
perl nFromAkw.pl -f InputDollarized.inp -p -z -m 1.16 outMinus.spectrum outPlus.spectrum
The mu is optional; if not given weights aren't computed.
Moreover, outMinus.spectrum can be outMinus1.spectrum,outMinus2.spectrum,...
and similarly for outPlus.spectrum. This is useful if you have multiple patches to build
a full omega, but not that patches of omega should NOT overlap for one sign (minus or plus). 
=cut

my $usage = "USAGE: $0 -f dollarizedInput [-m mu] [-p] [-z] minusSpectrum plusSpectrum\n";

my $templateInput;
my $isPeriodic;
my $zeroAtCenter = 0;
my $nonNegativeOnly = 0;
my $mu;

GetOptions('f=s' => \$templateInput,
           'm=f' => \$mu,
           'p' => \$isPeriodic,
		   'N' => \$nonNegativeOnly,
           'z' => \$zeroAtCenter) or die "$usage\n";

(defined($templateInput) and defined($isPeriodic)) or die "$usage\n";

my $sites;
my $eta;
my $geometryName;
my $geometryLeg = 1;
my $hptr = {"GeometryKind" => \$geometryName,
            "LadderLeg" => \$geometryLeg,
	    "TotalNumberOfSites" => \$sites,
            "CorrectionVectorEta" => \$eta};
OmegaUtils::getLabels($hptr,$templateInput);
die "$0 doesn't support LadderLeg=$geometryLeg\n" if ($geometryLeg > 2);

my ($fmString, $fpString) = @ARGV;
defined($fmString) or die "$usage\n";
defined($fpString) or print STDERR "$0: WARNING: Only one spectrum\n";

my (@filesMinus, @filesPlus);

getFiles(\@filesMinus, $fmString);
getFiles(\@filesPlus, $fpString) if (defined($fpString));

my $numberKs;
my %specMinus;
for (my $i = 0; $i < scalar(@filesMinus); ++$i) {
	readSpectrum(\%specMinus, \$numberKs, $filesMinus[$i]);
}

my %specPlus;
for (my $i = 0; $i < scalar(@filesPlus); ++$i) {
	readSpectrum(\%specPlus, \$numberKs, $filesPlus[$i]);
}

my %specFull;
addSpectrum(\%specFull, \%specMinus);
addSpectrum(\%specFull, \%specPlus) if (defined($fpString));
my $geometry = {"name" => $geometryName, "leg" => $geometryLeg};
OmegaUtils::printGnuplot(\%specFull, $geometry,  $isPeriodic, $zeroAtCenter, $nonNegativeOnly);
OmegaUtils::printOffsetPlots("offset", \%specFull, $geometry,  $isPeriodic, $zeroAtCenter);

exit(0) if (!defined($fpString));

my @nkx0;
my $norm = sumOverOmega(\@nkx0, \%specMinus, 0);
print "Norm=$norm\n";
if ($geometry->{"name"} eq "ladder") {
	my @nkxpi;
	$norm += sumOverOmega(\@nkxpi, \%specMinus, 1);
	printVsQ("outnkxpi.dat", \@nkxpi, $norm*$eta);
}

print "Norm=$norm\n";
printVsQ("outnkx0.dat", \@nkx0, $norm*$eta);

my $totalMy = ($geometry->{"name"} eq "ladder") ? 2 : 1;

my %fullVsOmega;
for (my $mp = 0; $mp < 2; ++$mp) { #mp = 0 is -, mp=1 is +
	my $ptr = ($mp == 0) ? \%specMinus : \%specPlus;
	for (my $my = 0; $my < $totalMy; ++$my) {
		my %h;
		sumOverKx(\%h, $ptr, $my);
		printVsOmega("nVsOmegaky$my"."Sector$mp.dat", \%h, $norm*$eta);
		addToFullOmega(\%fullVsOmega, \%h);
	}
}

printVsOmega("nVsOmega.dat", \%fullVsOmega, $norm*$eta);

if (defined($mu)) {
	sumWeight(\%fullVsOmega, $mu, $norm*$eta);
}

sub sumWeight
{
	my ($ptr, $mu, $scale) = @_;
	my ($below, $above) = (0, 0);
	my $max = 0;
	for my $omega (sort {$a <=> $b} keys %$ptr) {
		my $val = $ptr->{$omega};
		$max = $val if ($val > $max);
		#$val = 0 if ($val < 0);
		if ($omega < $mu) {
			$below += $val;
		} else {
			$above += $val;
		}
	}

	$max /= $scale;

	my $fout = "mu.dat";
	open(FOUT, ">", "$fout") or die "$0: Cannot write to $fout : $!\n";
	print FOUT "$mu 0\n";
	print FOUT "$mu $max\n";
	close(FOUT);
	print STDERR "File $fout written\n";

	my $factor = $below + $above;
	print STDERR "Factor= $factor, ";
	print STDERR "factor*eta/sites= ".$factor*$eta/$sites."\n";
	$factor = $sites/$factor;
	$below *= $factor;
	$above *= $factor;
	print STDERR "Below $mu : $below, above $mu: $above\n";
}

sub addToFullOmega
{
	my ($v, $ptr) = @_;
	for my $omega (sort {$a <=> $b} keys %$ptr) {
		my $val = $ptr->{$omega};
		$val = 0 if ($val < 0);
		if (!defined($v->{$omega})) {
			$v->{$omega} = $val;
		} else {
			$v->{$omega} += $val;
		}
	}
}

sub printVsOmega
{
	my ($fout, $ptr, $norm) = @_;
	open(FOUT, ">", "$fout") or die "$0: Cannot write to $fout : $!\n";
	for my $omega (sort {$a <=> $b} keys %$ptr) {
		my $val = $ptr->{$omega}/$norm;
		print FOUT "$omega $val\n";
	}

	close(FOUT);
	print STDERR "$0: File $fout has been written.\n";
}

sub sumOverKx
{
	my ($v, $ptr, $my) = @_;
	my ($factor, $fileIndices, $leg) = OmegaUtils::getGeometryDetails($geometry, $my);
	my $fileIndex = $my;

	for my $omega (sort {$a <=> $b} keys %$ptr) {
		my $aptr = $ptr->{$omega};
		my $nks = scalar(@$aptr) - 1;
		my $numberOfQs = int($factor*$nks);
		for (my $m2 = 0; $m2 < $numberOfQs; ++$m2) {
			#my $realPart = $aptr->[2*$m2+1+2*$fileIndex*$numberOfQs];
			my $imagPart = $aptr->[2*$m2+2+2*$fileIndex*$numberOfQs];
			if (defined($v->{$omega})) {
				$v->{$omega} += $imagPart;
			} else {
				$v->{$omega} = $imagPart;
			}
		}
	}
}

sub printVsQ
{
	my ($fout, $v, $norm) = @_;
	my $numberOfQs = scalar(@$v);
	my $centerShift = ($numberOfQs & 1) ? ($numberOfQs - 1)/2 : $numberOfQs/2;
	$centerShift = 0 unless ($zeroAtCenter);

	open(FOUT, ">", "$fout") or die "$0: Cannot write to $fout : $!\n";
	for (my $m2 = 0; $m2 < $numberOfQs; ++$m2) {
		my $m = $m2 - $centerShift;
		$m += $numberOfQs if ($m < 0);
		my $q = getQ($m2 - $centerShift, $numberOfQs, $isPeriodic);
		my $val = pi*$v->[$m]/$norm;
		print FOUT "$q $val\n";
	}

	close(FOUT);
	print STDERR "$0: File $fout has been written.\n";
}

sub sumOverOmega
{
	my ($v, $ptr, $my) = @_;
	my ($factor, $fileIndices, $leg) = OmegaUtils::getGeometryDetails($geometry, $my);

	my $fileIndex = $my;
	my $norm = 0;
	for my $omega (sort {$a <=> $b} keys %$ptr) { #no need to sort
		my $aptr = $ptr->{$omega};
		my $nks = scalar(@$aptr) - 1;
		my $numberOfQs = int($factor*$nks);
		for (my $m2 = 0; $m2 < $numberOfQs; ++$m2) {
			#my $realPart = $aptr->[2*$m2+1+2*$fileIndex*$numberOfQs];
			my $imagPart = $aptr->[2*$m2+2+2*$fileIndex*$numberOfQs];
			$norm += $imagPart;
			if (defined($v->[$m2])) {
				$v->[$m2] += $imagPart;
			} else {
				$v->[$m2] = $imagPart;
			}
		}
	}

	return $norm;
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

sub addSpectrum
{
	my ($v, $ptr) = @_;
	for my $omega (sort {$a <=> $b} keys %$ptr) { #no need to sort
		my $oldVal = $ptr->{$omega};
		my $aptr = $v->{$omega};
		if (defined($aptr)) {
			my $n = scalar(@$aptr);
			for (my $i = 1; $i < $n; ++$i) {
				$v->{$omega}->[$i] += $oldVal->[$i];
			}
		} else {
			my @temp = @$oldVal;
			$v->{$omega} = \@temp;
		}
	}
}


sub getQ
{
	my ($m, $n, $isPeriodic) = @_;
	return ($isPeriodic) ? 2.0*$pi*$m/$n : ($m + 1)*$pi/($n+1.0);
}

sub getFiles
{
	my ($fm, $string) = @_;
	my @temp = split(/,/, $string);
	@$fm = @temp;
}

