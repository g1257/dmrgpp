#!/usr/bin/perl -w

use strict;
use warnings;
use lib ".";
use Metts;

my ($beta,$label)=@ARGV;

my ($average, $total) = Metts::density($beta, $label, 1, *STDIN);
print STDERR "#Average= ".$average."\n";
print STDERR "#Total= $total\n";

