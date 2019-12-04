#!/usr/bin/perl

use strict;
use warnings;
use lib "../scripts";
use timeObservablesInSitu;


my ($site, $label, $tau)=@ARGV;
defined($tau) or die "USAGE: $0 site label tau < file\n";
my @matrix = timeObservablesInSitu::getMatrix($site, $label, *STDIN, $tau);
print "#site= $site\n";
print "#label=$label\n";
my $times = scalar(@matrix);

print STDERR "$0: Found $times times in STDIN\n";

for (my $t = 0; $t < $times; $t += 1) {
	my $time = $tau*$t;

	my $a = $matrix[$site];
	my $pair = $a->[$t];
	if (!defined($pair)) {
		next;
	}

	print "$time $site ".$pair->{"value"}." ".$pair->{"superdensity"}."\n";
}


