#!/usr/bin/perl

use strict;
use warnings;
use utf8;

my ($n, $m) = @ARGV;
defined($m) or die "Usage: $0 sites electronsUp\n";

my $total = combinatorial($n, $m);
print "TOTAL $total\n";
my $sum = 0;

# with contributions from Alberto N.
# N = number of sites
# Nup = Ndown = M
#numberOfDoubleOcc(D, N, M) = comb(N, D) * comb(N - D, M - D) * comb(N - M, M - D)
# comb(N, D) accounts for the number of ways of having doubly occupied sites on
# the lattice of N sites. Now, we have M-D ups and downs left to put somewhere.
#
#comb(N - D, M - D) accounts for the number of ways of putting the remaining M-D ups
# on the remaining N-D sites.
# Now, we need to fill the downs in the N-(M-D)-D=N-M spaces left, so
#
# comb(N-M, M - D) accounts for the number of ways of putting the
# remaining M-D downs on the remaining N-M sites left without ups.

for (my $d = 1; $d <= $m; ++$d) {
	my $r0 = factorialOf($n - $d);
	my $r1 = factorialOf($m - $d);
	my $r2 = factorialOf($n + $d - 2*$m);
	my $den = $r1*$r1*$r2;
	my $num = $d*$r0*combinatorial($n, $d);
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

