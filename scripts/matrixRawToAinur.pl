#!/usr/bin/perl

use strict;
use warnings;
use utf8;

my ($file, $label) = @ARGV;
defined($file) or die "USAGE: $0 filename [label]\n";
defined($label) or $label = "RAW_MATRIX";

open(FILE, "<", "$file") or die "$0: Cannot open $file : $!\n";
while (<FILE>) {
	last if (/^\Q$label/);

}

$_ = <FILE>;
chomp;
my @temp = split;
scalar(@temp) == 2 or die "$0: Expecting two integers, not $_\n";
my $n = $temp[0];
($n == $temp[1]) or die "$0: Expecting rows == cols, not $_\n";

print "$label=[\n";
for (my $row = 0; $row < $n; ++$row) {
	$_ = <FILE>;
	chomp;
	@temp = split;
	scalar(@temp) == $n or die "$0: Expecting row $row to contain $n columns\n";
	print "[".$temp[0];
	for (my $col = 1; $col < $n; ++$col) {
		print ", ".$temp[$col];
	}

	print "]";
	print "," if ($row + 1 < $n);
	print "\n";
}

close(FILE);
print "];\n";

