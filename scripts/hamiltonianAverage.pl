#!/usr/bin/perl -w
use strict;
use warnings;

my %val;
my %seenTime;
my $time = 0;

while(<STDIN>) {
	if (/Hamiltonian average at time=([^ ]+) for target=0 sector=[^ ]+ \<phi\(t\)\|H\|phi\(t\)\>=\(([^,]+),/) {
		$time = $1;
		if ($seenTime{$time}) {
			$val{$time}=0;
			$seenTime{$time}=0;
		}
		if (!defined($val{$time})) {
			$val{$time} = $2;
		} else {
			$val{$time} += $2;
		}
	} else {
		$seenTime{$time}=1;
	}
}

foreach my $key (sort {$a <=> $b} keys %val) {
	print "$key $val{$key}\n";
}

