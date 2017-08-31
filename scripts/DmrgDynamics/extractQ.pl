#!/usr/bin/perl

use strict;
use warnings;

my $q = shift @ARGV;
defined($q) or die "USAGE: $0 q file1 file2 ...\n";
my @files = @ARGV;

foreach my $file (@files) {
	procFile($file,$q);
}

sub procFile
{
	my ($file,$q) = @_;
	open(FILE, "<", $file) or die "$0: Cannot open file $file : $!\n";

	my $omega;
	while(<FILE>) {

		chomp;
		if (/^#omega=(.*$)/) {
			$omega = $1;
			next;
		}

		my @temp = split;
		next unless (scalar(@temp) == 3);
		next unless (abs($temp[0]-$q)<1e-3);
		die "$0: File $file line $_\n" if (!defined($omega));
		print "$omega $temp[1] $temp[2]\n";
		last;
	}

	close(FILE);
}

