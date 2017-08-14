#!/usr/bin/perl

use strict;
use warnings;
use Metts;

my ($beta,$betaLabel)=@ARGV;
my ($sum, $counter) = Metts::energy($beta, $betaLabel);
print STDERR "#Energy=$sum $counter\n";


