#!/usr/bin/perl -w

use strict;
use warnings;

my ($beta,$operator)=@ARGV;

defined($operator) or die "USAGE: $0 beta operator\n";

my $minSite = 1000;
my $maxSite = 0;
my @value;
my @value2;
my @counter;

my $label = "P0|".$operator;
while(<STDIN>) {
	if (/\Q$label/) {
		my @temp=split;
		my $site = $temp[0];
		my $val = $temp[1];
		my $t = $temp[2];
		next unless ($t==$beta);

		if (!defined($counter[$site])) {
			$counter[$site]=0;
		} else {
			$counter[$site]++;
		}
		$value[$site][$counter[$site]]=$val;
		$minSite = $site if ($site<$minSite);
		$maxSite = $site if ($site>$maxSite);
	}
}

my $total = $counter[1];
my ($denominator,$average)=(0,0);
defined($total) or die "$0: No data found\n";
for (my $i=0;$i<$total;$i++) {
	my $sum = 0;
	my $exitHere = 0;
	for (my $site=$minSite;$site<=$maxSite;$site++) {
		$_ = $value[$site][$i];
		if (!defined($_)) {
			$exitHere=1;
			last;
		}
		$sum += $_;
	}

	last if ($exitHere);


	print "$i $sum        ";
	$average += $sum;
	$denominator++;
	for (my $site=$minSite;$site<=$maxSite;$site++) {
		$_ = $value[$site][$i];
		$_ = 0 if (!defined($_));
		print "$_ ";
	}
	print "\n";
}
($denominator>0) or die "$0: No data found yet\n";
$average /= $denominator;
print STDERR "#Average= ".$average."\n";
print STDERR "#Total= $total\n";

