#!/usr/bin/perl

use strict;
use warnings;

package Ci;

sub getTests
{
	my $desc = "inputs/descriptions.txt";
	my ($tests) = @_;
	my $counter = 0;
	open(FILE,"$desc") or die "$0: Cannot open $desc : $!\n";

	while (<FILE>) {
		last if (/^\#TAGSTART/);
	}

	while (<FILE>) {
		last if (/^\#TAGEND/);
		next if (/^\#/);
		if (/^(\d+)\)/) {
			$tests->[$counter++] = $1;
		}
	}

	close(FILE);
}

sub isSu2
{
	my ($file,$n) = @_; 
	open(FILE, "$file") or return 0;
	my $su2 = 0;
	while (<FILE>) {
		chomp;
		if (/UseSu2Symmetry=1/) {
			$su2 = 1;
			last;
		}
	}

	close(FILE);
	return $su2;
}

1;

