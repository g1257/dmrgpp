#!/usr/bin/perl

use strict;
use warnings;
use utf8;

package EnergyAncillaInSitu;

sub main
{
	my ($beta, $betaLabel, $fin, $fout)=@_;
	defined($fout) or return;

	my $minSite = 1;
	my $minMeas = 0;

	my $counter  = 0;
	my $sum = 0;
	my $sum2 = 0;
	my $meas =0;
	my $site = -1;
	my $flag=0;

	while (<$fin>) {
		chomp;
		if (/sites=([^\+]+)\+([^\+]+)/) {
			my $site1 = $1;
			my $site2 = $2;
			$site = $site1;
			$site-- if ($site1>$site2);
			$meas++ if ($site==$minSite);
		}

		if (/Hamiltonian average at $betaLabel\=([^ ]+) for target\=0 /) {
			next unless ($1==$beta);
			next unless ($meas>$minMeas);
			if (/\<phi\(t\)\|H\|phi\(t\)>=/) {
				my @temp = split;
				my $tempSize = scalar(@temp);
				$tempSize > 1 or die "$0: Wrong line $_\n";
				my $energy = $temp[$tempSize-1];
				my ($re,$im) = complexCtor($energy);
				return "$0: $energy is not real\n" unless (abs($im)<1e-6);
				$sum2 += $re;
				$flag=1;
			}
		}

		if (/Hamiltonian average at $betaLabel\=([^ ]+) for target\=1 /) {
			next unless ($flag);
			$flag=0;
			$sum += $sum2;
			print $fout "$counter $sum2 $site\n";
			$sum2=0;
			$counter++;
		}
	}

	return "$0: counter==0\n" if ($counter == 0);

	$sum /= ($counter);
	print $fout "#Energy=$sum $counter\n";
	print $fout "#".$betaLabel."=$beta\n";
	return "";
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


