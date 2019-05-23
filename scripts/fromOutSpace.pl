#!/usr/bin/perl

use strict;
use warnings;
use Math::Trig;
use lib ".";
use OmegaUtils;
my ($outSpace, $templateInput, $isPeriodic, $zeroAtCenter) = @ARGV;
defined($templateInput) or die "USAGE: $0 out.space inputFile [isPeriodic] [zeroAtCenter]\n";
defined($isPeriodic) or $isPeriodic = 0;
defined($zeroAtCenter) or $zeroAtCenter = 0;

my $geometryName;
my $geometrySubName = "";
my $geometryLeg = 1;
my $centralSite;
my $hptr = {"GeometryKind" => \$geometryName,
            "GeometrySubKind" => \$geometrySubName,
            "LadderLeg" => \$geometryLeg,
            "TSPSites 1" => \$centralSite};

OmegaUtils::getLabels($hptr,$templateInput);

my $geometry = {"name" => $geometryName, "leg" => $geometryLeg, "subname" => $geometrySubName};

$hptr->{"centralSite"} = $centralSite;

my %h;
readSpaceValues(\%h, $outSpace);

my $outSpectrum = $outSpace;
$outSpectrum =~ s/\.space/\.spectrum/;
($outSpectrum ne $outSpace) or die "$0: $outSpectrum eq $outSpace ERROR FATAL\n";

open(FOUTSPECTRUM, ">", "$outSpectrum") or die "$0: Cannot write to $outSpectrum : $!\n";

foreach my $omega (sort keys %h) {
	print FOUTSPECTRUM "$omega ";
	my $spaceValues = $h{"$omega"};
	defined($spaceValues) or last;
	my @qValues;
	OmegaUtils::fourier(\@qValues,$spaceValues,$geometry,$hptr);
	my @array;
	OmegaUtils::writeFourier(\@array,\@qValues,$geometry);
	printSpectrum(\@array);
}

close(FOUTSPECTRUM);

print STDERR "$0: Wrote $outSpectrum\n";

sub readSpaceValues
{
	my ($h, $file) = @_;
	my $counter = 0;
	open(SPACEIN, "<", $file) or die "$0: Cannot open $file : $!\n";
	while (<SPACEIN>) {
		chomp;
		my ($omega, $n) = split;
		defined($n) or last;

		my @array;
		for (my $i = 0; $i < $n; ++$i) {
			$_ = <SPACEIN>;
			chomp;
			my ($i, $vv1, $vv2) = split;
			my @a = ($vv1, $vv2);
			$array[$i] = \@a;
		}

		$h->{"$omega"} = \@array;
		++$counter;
	}

	close(SPACEIN);
	print STDERR "$0: Read $counter omega values from $file\n";
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

