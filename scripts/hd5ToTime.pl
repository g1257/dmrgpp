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

	$_ = <PIPE>;
	print "$file $_";
	close(PIPE);
}

