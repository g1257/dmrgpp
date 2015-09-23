#!/usr/bin/perl

use strict;
use warnings;
use Math::Trig;
use OmegaUtils;

my ($omega0,$total,$omegaStep,$templateInput,$centralSite,$q) = @ARGV;
my $usage = "omegaBegin omegaTotal omegaStep dollarizedInput centralSite [q]\n";
$usage .= "\t q is the value of the momentum if positive or minus ";
$usage .= "the site if negative\n";
defined($centralSite) or die "USAGE: $0 $usage";


if ($omegaStep < 0) {
	my $beta = -$omegaStep;
	print STDERR "$0: Matsubara freq. assumed with beta= $beta\n";
	$omega0 = $omegaStep = 2.0*pi/$beta;
}

my $outSpectrum = "out.spectrum";
open(FOUTSPECTRUM,"> $outSpectrum") or die "$0: Cannot write to $outSpectrum : $!\n";
for (my $i = 0; $i < $total; ++$i) {

	my $omega = $omega0 + $omegaStep * $i;
	print FOUTSPECTRUM "$omega ";
	print STDERR "$0: About to proc for omega = $omega\n";

	if (defined($q)) {
		procThisOmega($i,$omega,$centralSite,$q);
	} else {
		my @array;
		procAllQs(\@array,$i,$omega,$centralSite);

		printSpectrum(\@array);
	}

	print STDERR "$0: Finished         omega = $omega\n";
}

close(FOUTSPECTRUM);
print STDERR "$0: Spectrum written to $outSpectrum\n";


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

	#my $cmd = "perl correctionVector.pl $omega < out$ind.txt > out$ind.space";
	#execThis($cmd);

	$inFile = "out$ind.space";
	$outFile = "out$ind.sq";

	my @spaceValues;
	readSpace(\@spaceValues,$inFile);
	my @qValues;
	my $geometry = getGeometry($templateInput);
	fourier(\@qValues,\@spaceValues,$geometry);
	printFourier($outFile,\@qValues,$geometry);
	#$cmd = "perl ft.pl $centralSite < out$ind.space > out$ind.sq";
	#execThis($cmd);
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

sub procThisOmega
{
	my ($ind,$omega,$centralSite,$q) = @_;
	procCommon($ind,$omega,$centralSite);

	if ($q < 0) {
		my $site = -$q;
		my $inFile = "out$ind.space";
		extractQ($inFile,$site);
		#$cmd = "perl extractQ.pl $site  out$ind.space";
	} else {
		my $inFile = "out$ind.sq";
		extractQ($inFile,$q);
		#$cmd = "perl extractQ.pl $q  out$ind.sq";
	}

	#execThis($cmd);
}

sub procAllQs
{
	my ($array,$ind,$omega,$centralSite) = @_;
	procCommon($ind,$omega,$centralSite);

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
		my $q = getQ($m,$n);
		my $ptr = $f->[$m];
		my @temp = @$ptr;
		my @f0 = @{$temp[0]};
		my @f1 = @{$temp[1]};
		my @sum = ($f0[0] - $f1[0],$f0[1] - $f1[1]);
		print FOUT "$q @sum\n";
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
	my $numberOfQs = $n;
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
	my $numberOfQs = int($n/2);
	for (my $m = 0; $m < $numberOfQs; ++$m) {
		my $q = getQ($m,$numberOfQs);
		my @f0 = fourierF0($v,$q);
		my @f1 = fourierF1($v,$q);
		my @sum = (\@f0,\@f1);
		$f->[$m] = \@sum;
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
		my $arg = $q*distanceLadder($i,$centralSite)*0.5;
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
	return pi*$m/($n+1.0);
	#return 2.0*pi*$m/$n;
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

sub extractQ
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


