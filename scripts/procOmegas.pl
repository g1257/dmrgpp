#!/usr/bin/perl

use strict;
use warnings;
use Math::Trig;
use OmegaUtils;
use Getopt::Long qw(:config no_ignore_case);

my $usage = "-b begin -t total -s step -f dollarizedInput -c ";
$usage .= "centralSite [-m mForQ] [-S site]\n";

my ($omega0,$total,$omegaStep,$templateInput,$centralSite,$site,$m);
my ($siteForSpectrum,$mForQ,$isPeriodic,$mMax,$wantsRealPart);

GetOptions('b=f' => \$omega0,
           't=i' => \$total,
		   's=f' => \$omegaStep,
		   'f=s' => \$templateInput,
		   'c=i' => \$centralSite,
		   'S:i' => \$siteForSpectrum,
		   'm:i' => \$mForQ,
		   'p' => \$isPeriodic,
		   'M:i' => \$mMax,
		   'r' => \$wantsRealPart) or die "$usage\n";

(
defined($omega0) and
defined($total) and
defined($omegaStep) and
defined($templateInput) and
defined($centralSite)
) or die "$0: USAGE: $usage\n";

if ($omegaStep < 0) {
	my $beta = -$omegaStep;
	print STDERR "$0: Matsubara freq. assumed with beta= $beta\n";
	$omega0 = $omegaStep = 2.0*pi/$beta;
}

my $geometry = getGeometry($templateInput);

my @omegas;
my $outSpectrum = "out.spectrum";
open(FOUTSPECTRUM,"> $outSpectrum") or die "$0: Cannot write to $outSpectrum : $!\n";
for (my $i = 0; $i < $total; ++$i) {

	my $omega = $omega0 + $omegaStep * $i;
	$omegas[$i] = $omega;
	print FOUTSPECTRUM "$omega ";
	print STDERR "$0: About to proc for omega = $omega\n";

	if (defined($mForQ)) {
		procThisOmegaKspace($i,$omega,$centralSite,$mForQ);
	} elsif (defined($siteForSpectrum)) {
		procThisOmegaSpace($i,$omega,$centralSite,$siteForSpectrum);
	} else {
		my @array;
		procAllQs(\@array,$i,$omega,$centralSite);

		printSpectrum(\@array);
	}

	print STDERR "$0: Finished         omega = $omega\n";
}

close(FOUTSPECTRUM);
print STDERR "$0: Spectrum written to $outSpectrum\n";
my $wantsRealOrImag = (defined($wantsRealPart)) ? "real" : "imag";
printSpectrumToColor($outSpectrum,$wantsRealOrImag,$geometry);
printGnuplot($outSpectrum,\@omegas,$geometry);

sub printGnuplot
{
	my ($inFile,$omegas,$geometry) = @_;

	open(FIN,"$inFile") or die "$0: Cannot open $inFile : $!\n";

	my @array;
	my $counter = 0;
	while (<FIN>) {
		my @temp = split;
		$array[$counter++] = \@temp;
	}

	close(FIN);

	my $numberOfOmegas = scalar(@$omegas);
	($counter == $numberOfOmegas) or
		die "$0: Found $counter omegas in $inFile but was expecting $numberOfOmegas\n";

	my $factor = 0;
	my @fileIndices=("");
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
		open(FOUT,"> $outFile") or die "$0: Cannot write to $outFile : $!\n";


		for (my $i = 0; $i < $numberOfOmegas; ++$i) {
			my $omega = $omegas->[$i];
			my $a = $array[$i];
			my $numberOfQs = int($factor*scalar(@$a));
			for (my $m = 0; $m < $numberOfQs; ++$m) {
				my $q = getQ($m,$numberOfQs);
				my $realPart = $a->[2*$m+$fileIndex*$numberOfQs];
				my $imagPart = $a->[2*$m+1];
				print FOUT "$q $omega $realPart $imagPart\n";
			}
		}

		close(FOUT);
		print "$0: Written $outFile\n";
	}
}

sub printSpectrumToColor
{
	my ($inFile,$what,$geometry) = @_;

	my @fileIndices=("");
	if ($geometry eq "chain") {
	} elsif ($geometry eq "ladder") {
		@fileIndices=(0,1);
	} else {
		die "$0: Unknown geometry $geometry\n";
	}

	foreach my $fileIndex (@fileIndices) {
		my $outSpectrum = "outSpectrum$fileIndex.color";
		my @colorData;
		my ($counter,$size) = spectrumToColor(\@colorData,
		                                      $inFile,
											  $what,
											  $geometry,
											  $fileIndex);

		open(FOUTSPECTRUM,"> $outSpectrum")
			or die "$0: Cannot write to $outSpectrum : $!\n";
		print FOUTSPECTRUM "$counter $size\n";

		my $rows = scalar(@colorData);
		for (my $i = 0; $i < $rows; ++$i) {
			my @thisRow = @{$colorData[$i]};
			my $cols = scalar(@thisRow);
			for (my $j = 0; $j < $cols; ++$j) {
				my $value = int($thisRow[$j]);
				print FOUTSPECTRUM $value." ";
			}

			print FOUTSPECTRUM "\n";
		}

		close(FOUTSPECTRUM);
		print STDERR "$0: Color spectrum written to $outSpectrum\n";
	}
}


sub printSpectrum
{
	my ($array) = @_;

	for (my $j = 0; $j < scalar(@$array); ++$j) {
		my $array2 = $array->[$j];
		my @array2 = @$array2;
		print FOUTSPECTRUM "$array2[1] $array2[2] ";
	}

	print FOUTSPECTRUM "\n";
}

sub procCommon
{
	my ($ind,$omega,$centralSite) = @_;
	my $n = OmegaUtils::getLabel($templateInput,"TotalNumberOfSites=");

	my $outFile = "out$ind.space";
	my $inFile = "out$ind.txt";
	my @values;
	my @values2;

	my $maxSite = correctionVectorRead(\@values,\@values2,$inFile);

	correctionVectorWrite($outFile,\@values,\@values2,$maxSite,$omega);

	$inFile = "out$ind.space";
	$outFile = "out$ind.sq";

	my @spaceValues;
	readSpace(\@spaceValues,$inFile);
	my @qValues;
	my $geometry = getGeometry($templateInput);
	fourier(\@qValues,\@spaceValues,$geometry);
	printFourier($outFile,\@qValues,$geometry);
}

sub correctionVectorRead
{
	my ($v1,$v2,$inFile) = @_;
	my $maxSite = 0;
	open(FIN,"$inFile") or die "$0: Cannot open $inFile : $!\n";
	my $status;
	while (<FIN>) {
		if (/P3/ and /PSI/) {
		        $status="p3";
		} elsif (/P2/ and /PSI/) {
		        $status="p2";
		} else {
		        next;
		}

		chomp;
		my @temp = split;
		die "$0: Line $_\n" unless (scalar(@temp)==5);

		my $site = $temp[0];
		$v1->[$site] = $temp[1] if ($status eq "p3");
		$v2->[$site] = $temp[1] if ($status eq "p2");
		$maxSite = $site if ($maxSite < $site);
		$status = "clear";
	}

	close(FIN);
	$maxSite++;

	print STDERR "$0: correctionVectorRead maxsite= $maxSite\n";
	return $maxSite;
}

sub correctionVectorWrite
{
	my ($outFile, $v1, $v2, $maxSite, $omega) = @_;

	open(FOUT,"> $outFile") or die "$0: Cannot write to $outFile : $!\n";

	print FOUT "#omega=$omega\n";
	for (my $i = 0; $i < $maxSite; ++$i) {
		my $vv1 = $v1->[$i];
		my $vv2 = $v2->[$i];
		defined($vv1) or die "$0: Undefined value for site = $i\n";
		print FOUT "$i $vv1 $vv2\n";
	}

	close(FOUT);
}

sub readSpace
{
	my ($space,$inFile) = @_;
	my $counter = 0;

	open(FIN,"$inFile") or die "$0: Cannot open $inFile : $!\n";
	while(<FIN>) {
		if (/^#/) {
		        next;
		}

		my @temp=split;
		next unless (scalar(@temp) == 3);
		my @temp2 = ($temp[1],$temp[2]);
		$space->[$temp[0]] = \@temp2;
		die "$0: Line $_\n" unless ($counter == $temp[0]);
		$counter++;
	}

	close(FIN);
	print STDERR "$0: Read $counter sites\n";
}

sub procThisOmegaKspace
{
	my ($ind,$omega,$centralSite,$mForQ) = @_;
	procCommon($ind,$omega,$centralSite);

	my $inFile = "out$ind.sq";
	extractValue($inFile,$mForQ);
}

sub procThisOmegaSpace
{
	my ($ind,$omega,$centralSite,$siteForSpectrum) = @_;
	procCommon($ind,$omega,$centralSite);

	my $inFile = "out$ind.space";
	extractValue($inFile,$siteForSpectrum);
}

sub procAllQs
{
	my ($array,$ind,$omega,$centralSite) = @_;
	procCommon($ind,$omega,$centralSite);
	readAllQs($array,$ind);
}

sub readAllQs
{
	my ($array,$ind) = @_;
	my $counter = 0;
	open(FILE,"out$ind.sq") or die "$0: Cannot open file : $!\n";
	while (<FILE>) {
		chomp;
		my @temp = split;
		my $n = scalar(@temp);
		($n == 3) or next;
		$array->[$counter++] = \@temp;
	}

	close(FILE);
}

sub execThis
{
	my ($cmd) = @_;
	print STDERR "$0: About to execute $cmd\n";
	system($cmd);
}

sub printFourier
{
	my ($outFile,$f,$geometry) = @_;
	if ($geometry eq "chain") {
		return printFourierChain($outFile,$f);
	}

	if ($geometry eq "ladder") {
		return printFourierLadder($outFile,$f);
	}

	die "$0: printFourier: undefined geometry $geometry\n";
}

sub printFourierChain
{
	my ($outFile,$f) = @_;

	open(FOUT,"> $outFile") or die "$0: Cannot write to $outFile : $!\n";

	my $n = scalar(@$f);
	for (my $m = 0; $m < $n; ++$m) {
		my $q = getQ($m,$n);
		my $ptr = $f->[$m];
		my @temp = @$ptr;
		print FOUT "$q $temp[0] $temp[1]\n";
	}

	close(FOUT);
}

sub printFourierLadder
{
	my ($outFile,$f) = @_;

	open(FOUT,"> $outFile") or die "$0: Cannot write to $outFile : $!\n";

	my $n = scalar(@$f);
	for (my $m = 0; $m < $n; ++$m) {
		my $ptr = $f->[$m];
		my @temp = @$ptr;
		print FOUT "$m @temp\n";
	}

	close(FOUT);
}

sub fourier
{
	my ($f,$v,$geometry) = @_;

	if ($geometry eq "chain") {
		return fourierChain($f,$v);
	}

	if ($geometry eq "ladder") {
		return fourierLadder($f,$v);
	}

	die "$0: ft: undefined geometry $geometry\n";
}

sub fourierChain
{
	my ($f,$v) = @_;
	my $n = scalar(@$v);
	my $numberOfQs = (defined($mMax)) ? $mMax : $n;
	for (my $m = 0; $m < $numberOfQs; ++$m) {
		my @sum = (0,0);
		my $q = getQ($m,$numberOfQs);
		for (my $i = 0; $i < $n; $i++) {
			my $ptr = $v->[$i];
			my @temp = @$ptr;
			my $arg = $q*($i-$centralSite);
			my $carg = cos($arg);
			my $sarg = sin($arg);
			$sum[0] += $temp[0]*$carg;
			$sum[1] += $temp[1]*$carg;
		}

		$f->[$m] = \@sum;
	}
}

sub fourierLadder
{
	my ($f,$v) = @_;
	my $n = scalar(@$v);
	my $numberOfQs = (defined($mMax)) ? $mMax : int(0.5*$n);
	for (my $m = 0; $m < $numberOfQs; ++$m) {
		my $q = getQ($m,$numberOfQs);
		my @f0 = fourierF0($v,$q);
		my @f1 = fourierF1($v,$q);
		for (my $x = 0; $x < 2; ++$x) {
			my $sign = 1-2*$x;
			my $realPart = $f0[0] + $sign*$f1[0];
			my $imagPart = $f0[1] + $sign*$f1[1];
			my @sum = ($realPart,$imagPart);
			$f->[$m+$numberOfQs*$x] = \@sum;
		}
	}
}

sub fourierF0
{
	my ($v,$q) = @_;
	my $n = scalar(@$v);
	my @sum;
	for (my $i = 0; $i < $n; $i+=2) {
		my $ptr = $v->[$i];
		my @temp = @$ptr;
		my $arg = $q*($i-$centralSite)*0.5;
		my $carg = cos($arg);
		$sum[0] += $temp[0]*$carg;
		$sum[1] += $temp[1]*$carg;
	}

	return @sum;
}

sub fourierF1
{
	my ($v,$q) = @_;
	my $n = scalar(@$v);
	my @sum;
	for (my $i = 1; $i < $n; $i+=2) {
		my $ptr = $v->[$i];
		my @temp = @$ptr;
		my $arg = $q*distanceLadder($i,$centralSite);
		my $carg = cos($arg);
		$sum[0] += $temp[0]*$carg;
		$sum[1] += $temp[1]*$carg;
	}

	return @sum;
}

sub distanceLadder
{
	my ($ind, $jnd) = @_;
	my $yind = ($ind & 1) ? 1 : 0;
	my $yjnd = ($jnd & 1) ? 1 : 0;
	my $xind = int($ind/2);
	my $xjnd = int($jnd/2);

	return ($xind - $xjnd);
}

sub getQ
{
	my ($m,$n) = @_;
	return ($isPeriodic) ? 2.0*pi*$m/$n : pi*$m/($n+1.0);
}

sub getGeometry
{
	my ($file) = @_;
	my $ret;
	my $label = "GeometryKind";
	open(FILE,$file) or die "$0: Cannot open $file : $!\n";
	while (<FILE>) {
		if (/$label=(.*$)/) {
			$ret = $1;
			last;
		}
	}

	close(FILE);

	defined($ret) or die "$0: Could not find $label in $file\n";
	return $ret;
}

sub extractValue
{
	my ($file,$q) = @_;
	open(FILE,$file) or die "$0: Cannot open file $file : $!\n";

	my $omega;
	while(<FILE>) {

		chomp;
		if (/^#omega=(.*$)/) {
			$omega = $1;
			next;
		}

		my @temp = split;
		next unless (scalar(@temp) > 1);
		next unless (abs($temp[0]-$q)<1e-3);
		die "$0: File $file line $_\n" if (!defined($omega));
		print "$omega $temp[1] ";
		print "$temp[2]" if (scalar(@temp) == 3);
		print "\n";
		last;
	}

	close(FILE);
}

sub spectrumToColor
{
	my ($data,$file,$realOrImag,$geometry,$qyIndex) = @_;
	my $counter = 0;
	my $size;
	my $finalSize;

	open(FIN,$file) or die "$0: Cannot open file $file : $!\n";
	while (<FIN>) {
		next if (/^#/);
		chomp;
		my @temp = split;
		my $n = scalar(@temp);
		next if ($n < 2);
		if (!defined($size)) {
			$size = $n;
		} else {
			($size == $n) or die "$0: Wrong line $_\n";
		}

		my @temp2 = getRealOrImagData(\@temp,$realOrImag,$geometry,$qyIndex);
		$finalSize = scalar(@temp2) if (!defined($finalSize));
		$data->[$counter++] = \@temp2;
	}

	close(FIN);

	print STDERR "$0: Read $counter lines size=$size for $realOrImag from $file\n";

	my ($min,$max) = minMaxData($data);
	print STDERR "$0: Data min = $min, max = $max\n";

	scaleData($data,$min,$max);
	return ($counter,$finalSize);
}

sub minMaxData
{
	my ($a) = @_;
	my $rows = scalar(@$a);
	my ($min,$max) = ($a->[0]->[1],$a->[0]->[1]);
	for (my $i = 0; $i < $rows; ++$i) {
		my @thisRow = @{$a->[$i]};
		my $cols = scalar(@thisRow);
		for (my $j = 0; $j < $cols; ++$j) {
			my $thisValue = $thisRow[$j];
			#next if ($thisValue<0);
			$min = $thisValue if ($min > $thisValue);
			$max = $thisValue if ($max < $thisValue);
		}
	}

	return ($min,$max);
}

sub scaleData
{
	my ($a,$min,$max) = @_;
	my $afactor = 255/($max-$min);
	my $bfactor = -$afactor*$min;
	my $rows = scalar(@$a);
	for (my $i = 0; $i < $rows; ++$i) {
		my @thisRow = @{$a->[$i]};
		my $cols = scalar(@thisRow);
		for (my $j = 0; $j < $cols; ++$j) {
			my $value = $a->[$i]->[$j];
			#if ($value < 0) {
			#	$a->[$i]->[$j] = 0;
			#	next;
			#}

			$a->[$i]->[$j] = $afactor*$value + $bfactor;
		}
	}
}

sub getRealOrImagData
{
	my ($d,$realOrImag,$geometry,$qyIndex) = @_;
	my @temp;
	my $n = scalar(@$d);
	my $start = 1;
	if ($geometry eq "ladder") {
		$n = int($n*0.5);
		$start += $qyIndex*$n;
		$n *= (1+$qyIndex);
	}

	my $j = 0;
	for (my $i = $start; $i < $n; ++$i) {
		next if ($realOrImag eq "imag" && ($i & 1));
		next if ($realOrImag eq "real" && !($i & 1));
		$temp[$j++] = $d->[$i];
	}

	return @temp;
}
