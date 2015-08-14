#!/usr/bin/env perl
#===============================================================================
#
#         FILE: timeEvolution3.pl
#
#        USAGE: ./timeEvolution3.pl
#
#  DESCRIPTION:
#
#      OPTIONS: ---
# REQUIREMENTS: ---
#         BUGS: ---
#        NOTES: ---
#       AUTHOR: G. A.
# ORGANIZATION:
#      VERSION: 1.0
#      CREATED: 08/14/2015
#     REVISION: ---
#===============================================================================

use strict;
use warnings;
use utf8;

while (<STDIN>) {
	chomp;
	my @temp = split;
	next if (/^#/);
	my $n = scalar(@temp);
	my $time = $temp[0];
	my $total = $n - 1;
	$total /= 2.0;
	my $index = 1;
	for (my $i = 0; $i < $total; ++$i) {
		my ($valr,$vali) = ($temp[$index],$temp[$index+1]);
		$index += 2;
		print "$time $i $valr $vali\n";
	}
}

