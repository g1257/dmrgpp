#!/usr/bin/perl

use strict;
use warnings;

package timeObservablesInSitu;

sub main
{
	my ($site, $label, $fin, $fout)=@_;
	defined($fout) or return;

	my $counter=0;
	print $fout "#site= $site\n";
	print $fout "#label=$label\n";

	while (<$fin>) {
		chomp;
		next unless /\Q$label/;
		if (/^${site} /) {
			my @temp = split;
			(scalar(@temp) >= 5) or next;
			my $value = procValue($temp[1]);
			my $time = $temp[2];
			my $superdensity = procValue($temp[4]);
			print $fout "$time  $value $superdensity\n";
		}
	}
}

sub procValue
{
	my ($t)=@_;
	$_=$t;
	s/\(//;
	s/,.*$//;
	return $_;
}

1;
