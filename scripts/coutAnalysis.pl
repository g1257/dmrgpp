#!/usr/bin/perl

=pod
USAGE is coutAnalysis.pl runForinput.cout [cutoff]
OUTPUT is 2 columns
Module TimeSpent
=cut

use warnings;
use strict;
use utf8;

my ($file, $cutoff) = @ARGV;
defined($file) or die "USAGE: $0 file\n";
defined($cutoff) or $cutoff = 0;

my %h;
my $totalTime = loadData(\%h, $file, $cutoff);

print STDERR "#cutoff=$cutoff\n";
print STDERR "#Elements=".scalar(keys %h)."\n";
printData(\%h);
print STDERR "#TotalForRun=$totalTime\n";

sub loadData
{
	my ($h, $file, $cutoff) = @_;
	open(FILE, "<", $file) or die "$0: Cannot open $file : $!\n";
	my $prevt = 0;
	my ($t0, $t1);
	while (<FILE>) {
		my $x = /^(.+) *\[(\d+)\]\:/;
		next unless ($x);
		my $name = $1;
		my $t = $2;
		$t0 = $t if (!defined($t0));
		$t1 = $t;
		my $dt = $t - $prevt;
		$prevt = $t;
		next if ($dt < $cutoff);
		if (defined($h->{$name})) {
			$h->{$name} += $dt;
		} else {
			$h->{$name} = $dt;
		}
	}

	close(FILE);
	return $t1 - $t0;
}


sub printData
{
	my ($hptr) = @_;
	my %h = %$hptr;
	my $tot = 0;
	foreach my $k (sort {$h{$b} <=> $h{$a}} keys %h) {
		my $t = $hptr->{$k};
		print "$k $t\n";
		$tot += $t;
	}

	print STDERR "--------------------\n";
	print STDERR "TotalForModules: $tot\n";
}

