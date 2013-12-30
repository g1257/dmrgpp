#!/usr/bin/perl

use strict;
use warnings;
use Utils;

my ($dmrgOrLanczos,$q,$root,$useReflection) = @ARGV;
defined($useReflection) or die "USAGE: $0 dmrgOrLanczos q root useReflection\n";
Utils::checkRange($dmrgOrLanczos,"Lanczos","Dmrg");

my $templateInput = "inputTemplate.inp";
my $b = 1;
my %params = Utils::loadParams("params.txt");

my $n = Utils::getLabel($templateInput,"TotalNumberOfSites=");
my @spectral;

for (my $site=0; $site<$n; $site++) {
	for (my $site2=0; $site2<$n; $site2++) {

		my ($siteMin,$siteMax) = findSubstitutes($site,$site2,$n,$useReflection);

		my $output = "$root${siteMin}_$siteMax";

		if ($b  && ($site2 >= $site)) {
			my ($psimagLite,$begin,$end,$step,$eps) = ($params{"PsimagLite"},$params{"OmegaBegin"},$params{"OmegaEnd"},
			                                           $params{"OmegaStep"},$params{"OmegaEps"});
			system("$psimagLite/drivers/continuedFractionCollection -f $output.comb -b $begin -e $end -s $step -d $eps > $output.cf");
		}


		if ($q == -2) {
			print STDERR "$0: Doing addDiagonal\n";
			addDiagonal("$output.cf",$site+$site2,$site,$site2,$n);
		} elsif ($q < 0) {
			addOptical("$output.cf",$site+$site2,$site,$site2,$n);
		} else {
			my @expq = findExp($site,$site2,$q);
			addThisCf("$output.cf",$site+$site2,\@expq);
		}

		print STDERR "$0: Done $site $site2";
		print STDERR " (using $siteMin $siteMax)" if ($siteMin != $site);
		print STDERR "\n";
	}
}

printCf("${root}Total$q.cf",\@spectral);

print "$0: Result is in ${root}Total$q.cf\n";

sub findSubstitutes
{
	my ($site,$site2,$n,$useReflection) = @_;

	my $siteMin = ($site < $site2) ? $site : $site2;
	my $siteMax = ($site < $site2) ? $site2 : $site;

	my ($s1,$s2) = ($siteMin,$siteMax);
	if ($useReflection and Utils::reflected($siteMin,$siteMax,$n)) {
		$s1 = Utils::findReflection($siteMax,$n);
		$s2 = Utils::findReflection($siteMin,$n);
	}

	return ($s1,$s2);
}

sub findExp
{
	my ($site,$site2,$q) = @_;
	my $momentum = ($q + 1) * acos(-1) / ($n + 1);
	#my $tmp =  ($site - $site2) * $momentum;
	#return (cos($tmp),sin($tmp));
	my $tmp = 2.0 * sin(($site+1)*$momentum) * sin(($site2+1)*$momentum) / ($n+1.0);
	return ($tmp,0);
}

sub acos { atan2( sqrt(1 - $_[0] * $_[0]), $_[0] ) }

sub addThisCf
{
	my ($file,$ind,$expq) = @_;
	open(FILE,"$file") or die "$0: Cannot open $file: $!\n";
	my $c = 0;
	while(<FILE>) {
		chomp;
		my @temp = split;
		next if (scalar(@temp) != 3);
		my $freq = $spectral[$c];
		my @freq;

		@freq = @$freq if ($ind > 0);
		checkFreq(\@freq,$ind,$temp[0]);

		my @value;
		$value[0] = $temp[1] * $expq->[0] + $temp[2] * $expq->[1];
		$value[1] = $temp[2] * $expq->[0] - $temp[1] * $expq->[1];

		$freq[1] += $value[0];
		$freq[2] += $value[1];

		$spectral[$c] = \@freq;
		$c++;
	}

	close(FILE);
}

sub addDiagonal
{
	my ($file,$ind,$site,$site2,$L) = @_;
	open(FILE,"$file") or die "$0: line $. Cannot open $file: $!\n";
	my $c = 0;
	return if ($site != $site2);

	my $diagFactor = ($site == $site2) ? 4 : 0.25;
	$diagFactor = 1 if ($dmrgOrLanczos eq "Lanczos");

	while(<FILE>) {
		chomp;
		my @temp = split;
		next if (scalar(@temp) != 3);
		my $freq = $spectral[$c];
		my @freq;

		@freq = @$freq if ($ind > 0);
		checkFreq(\@freq,$ind,$temp[0]);

		my @value;
		my $tmp = $diagFactor;
		$value[0] = $temp[1] * $tmp;
		$value[1] = $temp[2] * $tmp;

		$freq[1] += $value[0];
		$freq[2] += $value[1];

		$spectral[$c] = \@freq;
		$c++;
	}

	close(FILE);
}

sub addOptical
{
	my ($file,$ind,$site,$site2,$L) = @_;
	open(FILE,"$file") or die "$0: line $. Cannot open $file: $!\n";
	my $c = 0;
	my $diagFactor = ($site == $site2) ? 4 : 0.5;
	$diagFactor = 1 if ($dmrgOrLanczos eq "Lanczos");

	while(<FILE>) {
		chomp;
		my @temp = split;
		next if (scalar(@temp) != 3);
		my $freq = $spectral[$c];
		my @freq;

		@freq = @$freq if ($ind > 0);
		checkFreq(\@freq,$ind,$temp[0]);

		my @value;
		my $tmp = lMinusFactor($site,$L) * lMinusFactor($site2,$L) * $diagFactor;
		$value[0] = $temp[1] * $tmp;
		$value[1] = $temp[2] * $tmp;

		$freq[1] += $value[0];
		$freq[2] += $value[1];

		$spectral[$c] = \@freq;
		$c++;
	}

	close(FILE);
}

sub lMinusFactor
{
	my ($l,$L) = @_;
	return ($l + 1 - ($L+1)/2);
}


sub checkFreq
{
	my ($f,$ind,$omega) = @_;
	if ($ind == 0) {
		$f->[0] = $omega;
		$f->[1] = 0;
		$f->[1] = 0;
		return;
	}

	defined($f->[0]) or die "$0: Undefined freq for $omega\n";
	($f->[0] == $omega) or die "$0: Invalid frequency ".$f->[0]." should be $omega\n";
}

sub printCf
{
	my ($output,$f) = @_;
	my @freq = @$f;
	open(FOUT," > $output") or die "$0: Cannot open $output for writing: $!\n";
	for (my $i = 0; $i < scalar(@freq); $i++) {
		my $oneF = $freq[$i];
		my @oneFreq = @$oneF;
		print FOUT "$oneFreq[0] $oneFreq[1] $oneFreq[2]\n";
	}

	close(FOUT);
}

