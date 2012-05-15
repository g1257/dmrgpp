#!/usr/bin/perl -w
use strict;
use warnings;

my %val;
my %superd;
my %seenTime;
my $time = 0;

while(<STDIN>) {
	last if (/^ALL OPERATORS/);
}

while(<STDIN>) {
	if (/Hamiltonian average at time=([^ ]+) for target=0 sector=[^ ]+ \<phi\(t\)\|H\|phi\(t\)\>=\(([^,]+),[^\)]+\) +\<phi\(t\)\|phi\(t\)\>=\(([^,]+),/) {
		$time = $1;
		if ($seenTime{$time}) {
			$val{$time}=0;
			$superd{$time}=0;
			$seenTime{$time}=0;
		}
		if (!defined($val{$time})) {
			$val{$time} = $2;
		} else {
			$val{$time} += $2;
		}
		if (!defined($superd{$time})) {
			$superd{$time} = $3;
		} else {
			$superd{$time} += $3;
		}
	} else {
		$seenTime{$time}=1;
	}
}

print "#Time\t<H>\tSuperdensity\tQuotient\n";
foreach my $key (sort {$a <=> $b} keys %val) {
	my $quot = $val{$key}/$superd{$key};
	print "$key\t$val{$key}\t$superd{$key}\t$quot\n";
}

