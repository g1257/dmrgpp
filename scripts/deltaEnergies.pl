#!/usr/bin/perl

use strict;
use warnings;
use utf8;

my ($file) = @ARGV;
defined($file) or die "USAGE: $0 file\n";

my @energies;
my @ms;
loadEnergies(\@energies, \@ms, $file);
printDeltas(\@energies, \@ms);

sub loadEnergies
{
	my ($e, $m, $file) = @_;
	open(FILE, "<", $file) or die "$0: Cannot open $file : $!\n";
	while (<FILE>) {
		last if (/Infinite dmrg loop has been done/);
	}

	my $counter = 0;
	while (<FILE>) {
		chomp;
		if (/lowest eigenvalue= ([^ ]+)/) {
			$e->[$counter++] = $1;
		}

		if (/kept states=([^ ]+)/) {
			$m->[$counter] = $1;
		}
	}

	close(FILE);
}

sub printDeltas
{
	my ($e, $ms) = @_;
	my $n = scalar(@$e);
	return if ($n < 3);
	my $eprev = $e->[0];
	my $mprev = $ms->[0];
	for (my $i = 1; $i < $n; ++$i) {
		my $delta = $eprev - $e->[$i];
		$eprev = $e->[$i];
		my $m = $ms->[$i];
		defined($m) or $m = $mprev;
		print "$i $delta $m\n";
		$mprev = $m;
	}
}



