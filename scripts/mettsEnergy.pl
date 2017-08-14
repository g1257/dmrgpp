#!/usr/bin/perl

use strict;
use warnings;
use lib "../scripts";
use Metts;

my ($beta,$betaLabel)=@ARGV;
my ($sum, $counter) = Metts::energy($beta, $betaLabel,1,*STDIN);
print STDERR "#Energy=$sum $counter\n";


