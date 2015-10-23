#!/usr/bin/perl

use strict;
use warnings;

package OmegaUtils;

sub getLabels
{
	my ($hptr,$file) = @_;

	open(FILE,$file) or die "$0: Cannot open $file : $!\n";
	while (<FILE>) {
		chomp;
		foreach my $key (keys %$hptr) {
			if (/$key[= ]([^ ]+)/) {
				${$hptr->{$key}} = $1;
			}
		}
	}

	close(FILE);

	foreach my $key (keys %$hptr) {
		my $x = ${$hptr->{$key}};
		defined($x) or die "$0: Could not find $key in $file\n";
	}
}


1;

