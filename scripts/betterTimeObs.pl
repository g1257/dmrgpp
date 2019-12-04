#!/usr/bin/perl

use strict;
use warnings;
use utf8;
use timeObservablesInSitu;

my ($file, $label, $nsites, $tau) = @ARGV;
defined($nsites) or die "USAGE: $0 file label numberOfSites [tau]\n";
defined($tau) or $tau = 0.1;

my @matrix;
my $times = 0;
for (my $site = 0; $site < $nsites; ++$site) {

	open(FIN, "<", $file) or die "$0: Cannot open $file : $!\n";
	my @a = timeObservablesInSitu::getMatrix($site, $label, *FIN, $tau);
	$matrix[$site] = \@a;
	$_ = scalar(@a);
	$times = $_ if ($_ > $times);
	close(FIN);
}

print STDERR "$0: Found $times times in $file\n";

for (my $t = 0; $t < $times; $t += 1) {
	my $time = $tau*$t;
	for (my $site = 0; $site < $nsites; ++$site) {
		my $a = $matrix[$site];
		my $pair = $a->[$t];
		if (!defined($pair)) {
	#		print "$time $site -100 -100\n";
			next;
		}

		print "$time $site ".$pair->{"value"}." ".$pair->{"superdensity"}."\n";
	}
}


