#!/usr/bin/perl

use strict;
use warnings;
use utf8;

my ($n, $m) = @ARGV;
defined($m) or die "Usage: $0 sites electronsUp\n";

my $total = combinatorial($n, $m);
print "TOTAL $total\n";
my $sum = 0;
for (my $d = 1; $d <= $m; ++$d) {
	my $r0 = factorialOf($n - $d);
	my $r1 = factorialOf($m - $d);
	my $r2 = factorialOf($n + $d - 2*$m);
	my $den = $r1*$r1*$r2;
	my $num = $d*$d*$r0;
	my $val = $num / $den;
	$sum += $val;
}

$sum /= $total;
$sum /= $total;

print "$n $m $sum\n";

sub combinatorial
{
	my ($a, $b) = @_;
	my $sum = factorialOf($a)/factorialOf($b);
	return $sum/factorialOf($a - $b);
}

sub factorialOf
{
	my ($a) = @_;
	die "$0: factorialOf $a < 0\n" if ($a < 0);
	return 1 if ($a < 2);

	my $prod = 1;
	for (my $i = 1; $i <= $a; ++$i) {
		$prod *= $i;
	}

	return $prod;
}

