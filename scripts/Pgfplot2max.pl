#!/usr/bin/perl

use strict;
use warnings;
use utf8;

my ($file) = @ARGV;
defined($file) or die "USAGE: $0 filename\n";

my $pi = 3.14159;
my $points = 3;
my $kpoint = $pi;

my %max = maxFromFile($file);

my @disp = plotMax(\%max, 1);

my $slope = findSlope(\@disp, $kpoint, $points);

print STDERR "$0: Slope at $kpoint with $points points behind is $slope\n";

sub findSlope
{
	my ($a, $kpoint, $points) = @_;
	my $n = scalar(@$a);
	for (my $i = 0; $i < $n; ++$i) {
		my $ptr = $a->[$i];
		my ($k, $omega) = @$ptr;
		next if (abs($k - $kpoint) >= 1e-3);
		
		my $piIndex = $i;
		my $firstIndex = $i - $points;
		if ($firstIndex < 0) {
			die "$0: Points $points is too big\n";
		}
			
		my $ptr2 = $a->[$firstIndex];
		my ($k1, $omega1) = @$ptr2;
		my $slope = ($omega1 - $omega)/($k1 - $k);
		return $slope;
	}
}

sub plotMax
{
	my ($h, $print) = @_;
	my @disp;
	my $counter = 0;
	foreach my $k (sort {$a <=> $b} keys %$h) {
		my $ptr = $h->{$k};
		my ($omega, $value) = @$ptr;
		print "$k $omega\n" if ($print);
		$disp[$counter++] = [$k, $omega];
	}
	
	return @disp;
}

sub maxFromFile
{
	my ($file) = @_;
	my %max;
	open(FILE, "<", "$file") or die "$0: Cannot open $file : $!\n";
	while (<FILE>) {
		my @temp = split;
		next if (scalar(@temp) != 3);
		my $k = $temp[0];
		my $omega = $temp[1];
		my $value = $temp[2];
		
		my $ptr = $max{$k};
		my @newptr = ($omega, $value);
		if (defined($ptr)) {
			my ($oldOmega, $oldValue) = @$ptr;
			$max{$k} = \@newptr if ($value > $oldValue);
		} else {
			$max{$k} = \@newptr;
		}
	}
	
	close(FILE);
	
	return %max;
}
	
