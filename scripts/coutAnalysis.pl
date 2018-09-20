#!/usr/bin/perl

=pod
USAGE is coutAnalysis.pl runForinput.cout
OUTPUT is 2 columns
Module TimeSpent
=cut

use warnings;
use strict;
use utf8;

my ($file) = @ARGV;
defined($file) or die "USAGE: $0 file\n";

my %h;
my $totalTime = loadData(\%h, $file);

print STDERR "#Elements=".scalar(keys %h)."\n";
printData(\%h);
print STDERR "#TotalForRun=$totalTime\n";

sub loadData
{
	my ($h, $file) = @_;
	open(FILE, "<", $file) or die "$0: Cannot open $file : $!\n";

	my ($firstT, $lastT);
	while (<FILE>) {
		#LanczosSolver [7.714]: starting clock
		next unless (/^([^ ]+) \[([\d\.]+)\]\: ([^ ]+) clock/);
		my $name = $1;
		my $t = $2;
		my $what = $3;
		(($name) and ($t) and ($what)) or next;
		next if ($what ne "starting" and $what ne "stopping");

		$firstT = $t if (!defined($firstT));
		$lastT = $t;

		if (defined($h->{$name})) {
			my $ptr = $h->{$name};
			if (scalar(@$ptr) != 3) {
				print STDERR "$0: Error with $name\n";
				next;
			}

			my $value = $ptr->[1];
			my $dt = 0;
			if ($what eq "starting") {
				if ($value != 0) {
					print STDERR "$0: $name starting without stopping at line $.\n";
					next;
				}

				$value = $t;
			} elsif ($what eq "stopping") {
				if ($value > $t) {
					print STDERR "$0: $name has time skew\n";
					next;
				}

				$dt = $t - $value;
				$value = 0;
			} else {
				print STDERR "$0: expecting starting or stopping not $what\n";
				next;
			}
				
			my @temp = ($ptr->[0] + $dt, $value, $ptr->[2] + 1);
			$h->{$name} = \@temp;
		} else {
			if ($what ne "starting") {
				print STDERR "$0: expecting starting not $what\n";
				next;
			}

			my @temp = (0, $t, 1);
			$h->{$name} = \@temp
		}
	}

	close(FILE);
	return $lastT - $firstT;
}


sub printData
{
	my ($hptr) = @_;
	my %h = %$hptr;
	my $tot = 0;
	foreach my $k (sort {$h{$b}->[0] <=> $h{$a}->[0]} keys %h) {
		my $ptr = $hptr->{$k};
		if (scalar(@$ptr) != 3) {
			print STDERR "$0: Error with $k\n";
			next;
		}

		my $name = toFixedLength($k, 25, "after");
		my $sep = multiChar(" ",4);
		my $t = int($ptr->[0]*100)/100;
		my $time = toFixedLength($t, 6, "before");
		my $shouldBeZero = $ptr->[1];
		if ($shouldBeZero != 0) {
			print STDERR "$0: Error with $k, expecting 0, got $shouldBeZero\n";
			next;
		}

		my $numberOfTimes = toFixedLength($ptr->[2], 4, "before");
		print "$name$sep$time$sep$numberOfTimes\n";
		$tot += $t;
	}

	print STDERR "--------------------\n";
	print STDERR "TotalForObserved: $tot\n";
}

sub toFixedLength
{
	my ($what, $n, $spaceBeforeOrAfter) = @_;
	my $x = "$what";
	my $l = length($x);
	return $what if ($l >= $n);
	my $spaces = multiChar(" ", $n - $l);
	if ($spaceBeforeOrAfter eq "after") {
		return "$what$spaces";
	}
	
	return "$spaces$what";
}

sub multiChar
{
	my ($c, $n) = @_;
	my $sum = "";
	for (my $i = 0; $i < $n; ++$i) {
		$sum .= "$c";
	}

	return $sum;
}

