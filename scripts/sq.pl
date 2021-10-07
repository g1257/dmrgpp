#!/usr/bin/perl

use strict;
use warnings;
use utf8;
use Fourier;

my ($file, $label) = @ARGV;
defined($label) or die "USAGE: $0 filename label\n";

my @m = Fourier::loadMatrix($file, $label);

my @ft = Fourier::fourierTransform(\@m);

Fourier::printArray(\@ft);

