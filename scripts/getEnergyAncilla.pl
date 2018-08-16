#!/usr/bin/perl

use strict;
use warnings;
use utf8;
use lib "../scripts";
use EnergyAncillaInSitu;

my ($beta,$betaLabel)=@ARGV;
defined($beta) or die "USAGE: $0 beta [betaLabel] < file\n";

defined($betaLabel) or $betaLabel = "beta";

my $ret = EnergyAncillaInSitu::main($beta, $betaLabel, *STDIN, *STDOUT);
if ($ret ne "") {
	print STDERR "$ret\n";
}


