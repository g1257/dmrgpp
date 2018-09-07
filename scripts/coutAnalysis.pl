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
		my $x = /^(.+) *\[([\d\.]+)\]\:/;
		next unless ($x);
		my $name = $1;
		my $t = $2;
		$t0 = $t if (!defined($t0));
		$t1 = $t;
		my $dt = $t - $prevt;
		$prevt = $t;
		next if ($dt < $cutoff);

		if (defined($h->{$name})) {
			my $ptr = $h->{$name};
			if (scalar(@$ptr) != 2) {
				print STDERR "$0: Error with $name\n";
				last;
			}

			my @temp = ($ptr->[0] + $dt, $ptr->[1] + 1);
			$h->{$name} = \@temp;
		} else {
			my @temp = ($dt, 1);
			$h->{$name} = \@temp
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
	foreach my $k (sort {$h{$b}->[0] <=> $h{$a}->[0]} keys %h) {
		my $ptr = $hptr->{$k};
		if (scalar(@$ptr) != 2) {
			print STDERR "$0: Error with $k\n";
			last;
		}

		my $t = int($ptr->[0]);
		my $numberOfTimes = $ptr->[1];
		print "$k\t $t\t $numberOfTimes\n";
		$tot += $t;
	}

	print STDERR "--------------------\n";
	print STDERR "TotalForModules: $tot\n";
}

