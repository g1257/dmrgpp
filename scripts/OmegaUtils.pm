#!/usr/bin/perl

use strict;
use warnings;

package OmegaUtils;

sub getLabel
{
	my ($file,$label)=@_;
	open(FILE,"$file") or die "$0: Cannot open $file: $!\n";
	my $value;
	while (<FILE>) {
	        chomp;
	        if (/$label(.*$)/) {
	                $value=$1;
	                last;
	        }
	}

	close(FILE);

	defined($value) or die "$0: Could not find $label in $file\n";

	return $value;
}


1;
