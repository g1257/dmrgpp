#!/usr/bin/perl

use strict;
use warnings;
use Math::Trig;
use lib ".";
use OmegaUtils;
my ($outSpectrum, $templateInput, $isPeriodic, $zeroAtCenter) = @ARGV;
defined($templateInput) or die "USAGE: $0 out.spectrum inputFile [isPeriodic] [zeroAtCenter]\n";
defined($isPeriodic) or $isPeriodic = 0;
defined($zeroAtCenter) or $zeroAtCenter = 0;

my $geometryName;
my $geometrySubName = "";
my $geometryLeg = 1;
my $hptr = {"GeometryKind" => \$geometryName,
            "GeometrySubKind" => \$geometrySubName,
            "LadderLeg" => \$geometryLeg};

OmegaUtils::getLabels($hptr,$templateInput);

my $geometry = {"name" => $geometryName, "leg" => $geometryLeg, "subname" => $geometrySubName};

OmegaUtils::printGnuplot($outSpectrum, $geometry, $isPeriodic, $zeroAtCenter);
