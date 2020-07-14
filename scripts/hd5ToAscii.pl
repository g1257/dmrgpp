#!/usr/bin/perl

use strict;
use warnings;
use utf8;

my ($file) = @ARGV;
defined($file) or die "USAGE: $0 filename\n";

open(FILE, "<", $file) or die "$0: Cannot open $file : $!\n";

my $lastTitle;
my @buffer;
my %h;
my $lengthOfData;
my $sites;

while (<FILE>) {
	chomp;
	if (/\(0\)\: (.+$)/) {
		$sites = $1;
		last;
	}
}

while (<FILE>) {
	chomp;
	next if (/DATASPACE /);
	if (/DATASET (.+$)/) {
		my $n = scalar(@buffer);
		if ($n > 0 and !defined($lastTitle)) {
			die "$0: Saw (\\d+): before DATASET\n";
		}

		if ($n > 0) {
			if (defined($lengthOfData) and $lengthOfData != $n) {
				die "$0: Sets of different sizes for $lastTitle $lengthOfData != $n\n";
			}

			$lengthOfData = $n;

			my @copy = @buffer;
			$h{"$lastTitle"} = \@copy;
			@buffer = ();
		}

		$lastTitle=$1;
		$lastTitle=~s/\"//g;
		$lastTitle=~s/[\{\}]//g;
		$lastTitle=~s/ //g;
		next;
	}

	if (/\(\d+\)\:(.+$)/) {
		my $tmp = $1;
		$tmp =~ s/, *$//;
		my @temp = split/,/, $tmp;
		push @buffer, @temp;
		next;
	}
}

close(FILE);

my $n = scalar(@buffer);
if ($n > 0) {
	if (defined($lengthOfData) and $lengthOfData != $n) {
		die "$0: Sets of different sizes for $lastTitle $lengthOfData != $n\n";
	}

	my @copy = @buffer;
	$h{"$lastTitle"} = \@copy;
	@buffer = ();
}

printHash(\%h, $lengthOfData);

sub printHash
{
	my ($h, $lengthOfData) = @_;
	my $i = 0;
	for (my $l = 0; $l < $lengthOfData; ++$l) {
		print "$l";
		foreach my $key (sort {$a <=> $b} keys %$h) {
			my $ptr = $h->{"$key"};
			my $m = scalar(@$ptr);
			print " ".$ptr->[$l];
		}

		print "\n";
	}
}

