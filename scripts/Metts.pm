#!/usr/bin/perl

use strict;
use warnings;

package Metts;

sub energy
{
	my ($beta,$betaLabel,$option,$fin) = @_;
	defined($fin) or die "Metts::energy beta label option fin\n";

	my $minSite = 1;
	my $minMeas = 0;

	my $counter  = 0;
	my $sum = 0;
	my $sum2 = 0;
	my $meas =0;
	my $site = -1;
	my $flag=0;
	while (<$fin>) {
		if (/sites=([^\+])\+([^\+])/) {
			my $site1 = $1;
			my $site2 = $2;
			$site = $site1;
			$site-- if ($site1>$site2);
			$meas++ if ($site==$minSite);
		}

		if (/Hamiltonian average at $betaLabel\=([^ ]+) for target\=0 /) {
			next unless ($1==$beta);
			next unless ($meas>$minMeas);
			next unless ($site==$minSite);
		 	if (/\<phi\(t\)\|H\|phi\(t\)>=([^ ]+) /) {
				my $energy = $1;
				my ($re,$im) = complexCtor($energy);
				die "$0: $energy is not real\n" unless (abs($im)<1e-6);
				$sum2 += $re;
				$flag=1;
			}
		}

		if (/Hamiltonian average at $betaLabel\=([^ ]+) for target\=1 /) {
			next unless ($flag);
			$flag=0;
			$sum += $sum2;
			print "$counter $sum2 $site\n" if ($option);
			$sum2=0;
			$counter++;
		}
	}

	die "$0: counter==0\n" if ($counter==0);
	$sum /= ($counter);
	return ($sum, $counter);
}

sub density
{
	my ($beta,$label,$option,$fin)=@_;

	defined($fin) or die "Metts::density: beta label option fin\n";

	#print STDERR "$beta $label $option\n";
	my $minSite = 1000;
	my $maxSite = 0;
	my @value;
	my @value2;
	my @counter;

	while(<$fin>) {
		next unless (/\Q$label/);
		my @temp=split;
		(scalar(@temp) > 2) or next;

		my $site = $temp[0];
		my $val = $temp[1];
		my $t = $temp[2];
		next unless ($t == $beta);

		if (!defined($counter[$site])) {
			$counter[$site]=0;
		} else {
			$counter[$site]++;
		}

		$value[$site][$counter[$site]]=$val;
		$minSite = $site if ($site<$minSite);
		$maxSite = $site if ($site>$maxSite);
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

		print "$i $sum        " if ($option);
		$average += $sum;
		$denominator++;
		for (my $site=$minSite;$site<=$maxSite;$site++) {
			$_ = $value[$site][$i];
			$_ = 0 if (!defined($_));
			print "$_ "  if ($option);
		}

		print "\n"  if ($option);
	}

	($denominator>0) or die "$0: No data found yet\n";
	$average /= $denominator;
	return ($average, $total);
}


sub complexCtor
{
	my ($x) = @_;
	if ($x=~/\(([^,]+),([^\)]+)\)/) {
		return ($1,$2);
	}

	die "$0: $x is not numeric\n" unless ($x=~/^[0-9e\-\+\.]+$/);
	return ($x,0);
}

1;
