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
defined($cutoff) or $cutoff = 2;

my %h;
loadData(\%h, $file, $cutoff);

print STDERR "#Elements=".scalar(keys %h)."\n";
printData(\%h);

sub loadData
{
	my ($h, $file, $cutoff) = @_;
	open(FILE, "<", $file) or die "$0: Cannot open $file : $!\n";
	my $prevt = 0;
	while (<FILE>) {
		my $x = /^(.+) *\[(\d+)\]\:/;
		next unless ($x);
		my $name = $1;
		my $t = $2;
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
}


sub printData
{
	my ($hptr) = @_;
	my %h = %$hptr;
	foreach my $k (sort {$h{$b} <=> $h{$a}} keys %h) {
		my $t = $hptr->{$k};
		print "$k $t\n";
	}
}

