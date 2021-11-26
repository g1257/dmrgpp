#!/usr/bin/perl

use strict;
use warnings;
use Math::Trig;
use lib ".";
use OmegaUtils;
my ($outSpectrum, $templateInput, $isPeriodic, $zeroAtCenter, $nonNegativeOnly) = @ARGV;
defined($templateInput) or die "USAGE: $0 out.spectrum inputFile [isPeriodic] [zeroAtCenter] [nonNegativeOnly]\n";
defined($isPeriodic) or $isPeriodic = 0;
defined($zeroAtCenter) or $zeroAtCenter = 0;

my $isAinur = OmegaUtils::isAinur($templateInput);
my $geometryName;
my $geometrySubName = "";
my $geometryLeg = 1;
my $hptr = {"GeometryKind" => \$geometryName,
            "GeometrySubKind" => \$geometrySubName,
            "LadderLeg" => \$geometryLeg};

OmegaUtils::getLabels($hptr,$templateInput);

if ($isAinur) {
	$geometryName =~ s/\"//g;
	$geometryName =~ s/ *; *$//;
}


my $geometry = {"name" => $geometryName, "leg" => $geometryLeg, "subname" => $geometrySubName};

OmegaUtils::printGnuplot($outSpectrum, $geometry, $isPeriodic, $zeroAtCenter, $nonNegativeOnly);
