#!/usr/bin/perl

use strict;
use warnings;
use utf8;

my $cmd = "h5dump -d /Def/FinalPsi/TimeSerializer/Time ";

my @files = (@ARGV);

foreach my $file (@files) {
	my $cmd2 = "$cmd \"$file\"";
	my $ret = open(PIPE, "$cmd2 |");
	if (!$ret) {
		print STDERR "$0: Could not open pipe $cmd2\n";
		next;
	}

	my $time;
	while (<PIPE>) {
		chomp;
		if (/\(0\)\: (.*$)/) {
			$time = $1;
			last;
		}
	}

	if (!$time) {
		print STDERR "$0: Could not find Time in $file\n";
		next;
	}

	print "$file $time\n";
	close(PIPE);
}

