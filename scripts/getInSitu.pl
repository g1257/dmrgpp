#!/usr/bin/perl

use strict;
use warnings;
use utf8;

my ($file, $label, $lx, $ly) = @ARGV;
defined($lx) or die "USAGE: $0 filename label lx [ly]\n";
defined($ly) or $ly = 1;

my @value;
open(FILE, "<", $file) or die "$0: Cannot open $file : $!\n";

while (<FILE>) {
	if (/\Q$label/) {
	#63 0.899850 0.000000 <gs|n|gs> 1.000000
		my @temp = split;
		next if (scalar(@temp) != 5);
		my $obs = $temp[3];
		next unless ($obs eq $label);
		$value[$temp[0]] = $temp[1];
	}
}

close(FILE);

my $total = 0;
for (my $y = 0; $y < $ly; ++$y) {
	for (my $x = 0; $x < $lx; ++$x) {
		my $ind = $y + $x*$ly;
		my $tmp = $value[$ind];
		defined($tmp) or $tmp = "X";
		$total += $tmp unless ($tmp eq  "X");
		print "$tmp ";
	}

	print "\n";
}

print "#TOTAL=$total\n";

