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
tbw
=cut

my $usage = "USAGE: $0 -f dollarizedInput -m mu -k kx [-p] [-z] out.gnuplot\n";

my $templateInput;
my $isPeriodic;
my $zeroAtCenter = 0;
my $mu;
my $kx;

GetOptions('f=s' => \$templateInput,
           'm=f' => \$mu,
           'p' => \$isPeriodic,
           'z' => \$zeroAtCenter,
           'k=f' => \$kx) or die "$usage\n";

my $x = (defined($kx) & defined($mu));
($x and defined($templateInput) and defined($isPeriodic)) or die "$usage\n";

my $geometry;
my $sites;
my $eta;
my $hptr = {"GeometryKind" => \$geometry,
	    "TotalNumberOfSites" => \$sites,
            "CorrectionVectorEta" => \$eta};
OmegaUtils::getLabels($hptr,$templateInput);
my ($file) = @ARGV;

defined($file) or die "$usage\n";
my ($gamma, $delta, $anorm);
getMinParams(\$gamma, \$delta, \$anorm, $file, $mu, $kx);
$x = (defined($gamma) and defined($delta) and defined($anorm));
$x or die "$0: Fit failed\n";
print STDERR "gamma=$gamma delta=$delta anorm=$anorm\n";

my %spec;
my $numberOfKs;
my @omegas;
findOmegas(\@omegas, \$numberOfKs, $file);
createSpectrum(\%spec, \@omegas, $numberOfKs);
OmegaUtils::printOffsetPlots("fit", \%spec, $geometry,  $isPeriodic, $zeroAtCenter);

sub findOmegas
{
	my ($omegas, $numberOfKs, $file) = @_;
	open(FILE, "<", $file) or die "$0: Cannot open $file : $!\n";
	my $counter = 0;
	my $x = 0;
	my $omegaPrev;
	while (<FILE>) {
		my @temp = split;
		if (defined($omegaPrev) && $omegaPrev != $temp[1]) {
			$$numberOfKs = $counter;
			$counter = 0;
			$omegas[$x++] = $omegaPrev;
		}

		$omegaPrev = $temp[1];
		++$counter;
	}

	close(FILE);
	$omegas[$x++] = $omegaPrev;
	print STDERR "$0: Found ".scalar(@omegas)." omegas\n";
	print STDERR "$0: Found numberOfKs= ".$$numberOfKs."\n";
}

sub createSpectrum
{
	my ($ptr, $omegas, $numberOfKs) = @_;
	foreach my $omega (@$omegas) {
		my @temp = generateSpectrum($omega, $numberOfKs);
		$ptr->{$omega} = \@temp;
	}
}

sub generateSpectrum
{
	my ($omega, $numberOfQs) = @_;

	my $omega2 = $omega*$omega;
	my $gaom = $gamma*$omega;
	my @result;
	$result[0] = $omega;
	for (my $ky = 0; $ky < 2; ++$ky) {
		for (my $m2 = 0; $m2 < $numberOfQs; ++$m2) {
			my $q = OmegaUtils::getQ($m2, $numberOfQs, $isPeriodic);
			my $ek = dispersionBar($q, $ky*$pi);
			my $num = ($omega + $ek)*2.0*$gaom*$anorm/$pi;
			my $phi2 = $ek**2 + $gamma**2 + $delta**2;
			my $den = ($omega2 - $phi2)**2 + 4*$gaom*$gaom;
			my $imagPart = $num/$den;
			$result[2*$m2+1+2*$ky*$numberOfQs] = 0; # real part
			$result[2*$m2+2+2*$ky*$numberOfQs] = $imagPart;     # imag part
		}
	}

	return @result;
}

sub dispersionBar
{
	my ($kx, $ky) = @_;
	return -2*cos($kx) - cos($ky) - $mu;
}

sub getMinParams
{
	my ($gamma, $delta, $anorm, $file, $mu, $kx) = @_;
	my $ky = 0;
	my $cmd = "../../PsimagLite/drivers/fit $file $mu $kx $ky";
	open(PIPE, "$cmd |") or die "$0: Cannot open pipe : $!\n";
	while (<PIPE>) {
		next if (/^#/);
		if (/gamma=(.+$)/) {
			$$gamma = $1;
			next;
		}

		if (/delta=(.+$)/) {
			$$delta = $1;
			next;
		}

		if (/anorm=(.+$)/) {
			$$anorm = $1;
			next;
		}
	}

	close(PIPE);
}



