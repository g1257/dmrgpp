#!/usr/bin/perl

use strict;
use warnings;

my ($beta,$label) = @ARGV;
defined($label) or die "USAGE: $0 beta label < observeOutput\n";

my $label2 = "site $label time";
my $flag = 0;
my @value;
my $maxSite = 0;
my $counter = -1;
while (<STDIN>) {
	if (!$flag && /\Q$label2/) {
		$flag = 1;
		$counter++;
		next;
	}

	next if (!$flag);

	if (/^DOFS=/) {
		$flag = 0;
		next;
	}

	my @temp = split;

	if (scalar(@temp) != 3) {
		next;
	}

	if (!isNumeric($temp[2]) || $temp[2] != $beta) {
		next;
	}

	$value[$counter][$temp[0]]=$temp[1];
	$maxSite = $temp[0] if ($temp[0] > $maxSite);
}

print STDERR "$0: Found $counter vectors, maxSite=$maxSite\n";

for (my $i = 0; $i < $counter; ++$i) {
	print "$i ";
	my $undef = 0;
	my $sum = 0;
	for (my $site = 0; $site <= $maxSite; ++$site) {
		my $x = $value[$i][$site];
		if (!defined($x)) {
			$x = "UNDEF";
			$undef = 1;
		} else {
			$sum += $x;
		}

		print "$x ";
	}

	$sum = "UNDEF" if ($undef);
	print "$sum \n";
}

sub isNumeric
{
	my ($x) = @_;
	return 1 if ($x =~ /^[\d\.e\-+]+$/);

	return 0;
}

