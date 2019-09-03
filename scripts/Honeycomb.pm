#!/usr/bin/perl

use strict;
use warnings;
use utf8;
use Math::Trig;

package Honeycomb;

my $pi = Math::Trig::pi;

sub setSiteCoordinates
{
	my ($tindx, $tindy, $ptsx, $ptsy, $ncmatrix, $lx, $ly, $isHash) = @_;
	my $isArmchairX = $isHash->{"isArmchairX"};

	if ($isArmchairX){
		honeycombArmchairSetSiteCoordinates($tindx, $tindy, $ptsx, $ptsy, $ncmatrix, $lx, $ly, $isHash);
	} else {
		honeycombZigzagSetSiteCoordinates($tindx, $tindy, $ptsx, $ptsy, $ncmatrix, $lx, $ly, $isHash);
	}
}

sub honeyGetQ
{
	my ($m1, $m2, $lx, $ly, $type) = @_;
	my ($b1x, $b1y) = (2*$pi/(3*$lx), 0);
	my ($b2x, $b2y) = (0, 2*$pi/(sqrt(3)*$ly));
	if ($type eq "zigzag") {
		($b1x, $b1y) = (2*$pi*sqrt(3)/(3*$lx), -1);
		($b2x, $b2y) = (0, 4*$pi/(3*$ly));
	}

	return ($m1 * $b1x + $m2 * $b2x, $m1 * $b1y + $m2 * $b2y);
}

sub honeySpace
{
	my ($tindx, $tindy, $n, $hptr, $type) = @_;
	my %isHash = %$hptr;
	$isHash{"isArmchairX"} = ($type =~ /armchair/i);
	$isHash{"isLeft"} = ($hptr->{"#options"} =~ /left/i);
	my (@ptsx, @ptsy, $ncmatrix);
	my $lx = $hptr->{"#Lx"};
	my $ly = $hptr->{"#Ly"};
	setSiteCoordinates($tindx, $tindy, \@ptsx, \@ptsy, $ncmatrix, $lx, $ly, \%isHash);
}

sub fillQvalues
{
	my ($hptr, $type) = @_;

	my $lx = $hptr->{"#Lx"};
	my $ly = $hptr->{"#Ly"};
	my ($M1, $M2) = (2*$lx, $ly);
	my @array;	
	for (my $m1 = 0; $m1 < $M1; ++$m1) {
		for (my $m2 = 0; $m2 < $M2; ++$m2) {
			my ($qx, $qy) = honeyGetQ($m1, $m2, $lx, $ly, $type);
			# choose your path
			next if ($qy != 0);
			my $m = $m1 + $m2*$M1;
			my $temp = { "m" => $m, "q" => $qx};
			push @array, $temp;
		}
	}

	return @array;
}

sub honeycombZigzagSetSiteCoordinates
{
	my ($tindx, $tindy, $ptsx, $ptsy, $ncmatrix, $lx, $ly, $isHash) = @_;
	my $isLeft = $isHash->{"isLeft"};
	my ($scalex, $scaley) = (2.0, 4.0/sqrt(3.0));

	### right is the default
	### x-scale is scaled by 2 and y is scaled by 4/sqrt(3)
	my @t1 = (2, 0); #(1.0*$scalex, 0*$scaley); ### [2,0]
	my @t2 = (1, 2); #(0.5*$scalex,sqrt(3.0)/2.0 *$scaley); ### [1,2]

	if ($isLeft) {
		### x-scale is scaled by 2 and y is scaled by 4/sqrt(3)
		@t1 = (2, 0); #(1.0*$scalex, 0*$scaley); ### [2,0]
		@t2 = (1, 2); #(-0.5*$scalex, sqrt(3.0)/2.0 *$scaley); ### [1,2]
	}

	my $ncrows = $ncmatrix->{"rows"};
	my $counter = 0;

	for (my $i = 0; $i < $lx; ++$i) {
		my $xv = $i*$t1[0];
		my $txv = $i;
		my $x0 = $xv;
		my $y0 = 0;
		my $tx0 = $txv;
		my $ty0 = 0;
		my $sign = ($isLeft) ? -1 : 1;
		my $x1 = $xv + $sign;
		my $y1 = 1.0;
		my $tx1 = $txv + 0.5*$sign;
		my $ty1 = 0.5/sqrt(3.0);

		for (my $j = 0; $j < $ly; ++$j) {
			my $cxa = int($x0 + $j*$t2[0]);
			my $cxb = int($x1 + $j*$t2[0]);
			my $cya = int($y0 + $j*$t2[1]);
			my $cyb = int($y1 + $j*$t2[1]);

			my $txa = $tx0 + $j*$t2[0]/$scalex;
			my $txb = $tx1 + $j*$t2[0]/$scalex;
			my $tya = $ty0 + $j*$t2[1]/$scaley;
			my $tyb = $ty1 + $j*$t2[1]/$scaley;

			$tindx->[$counter] = $txa;
			$tindy->[$counter] = $tya;
			#print x0, y0, cxa, cya
			push @$ptsx, $cxa;
			push @$ptsy, $cya;
			$ncmatrix->{"data"}->[$cxa + $cya*$ncrows] = $counter++;

			$tindx->[$counter] = $txb;
			$tindy->[$counter] = $tyb;
			push @$ptsx, $cxb;
			push @$ptsy, $cyb;
			$ncmatrix->{"data"}->[$cxb + $cyb*$ncrows] = $counter++;
		}
	}
}


sub honeycombArmchairSetSiteCoordinates
{
	my ($tindx, $tindy, $ptsx, $ptsy, $ncmatrix, $lx, $ly, $isHash) = @_;

	my $isLeft = $isHash->{"isLeft"};
	my ($scalex, $scaley) = (1.0, 2.0/sqrt(3.0));

	### right is the default
	### x-scale is scaled by 2 and y is scaled by 4/sqrt(3)
	my @t1 = (3, 0); #(1.0*$scalex, 0*$scaley); ### [2,0]
	my @t2 = (0, 2); #(0.5*$scalex,sqrt(3.0)/2.0 *$scaley); ### [1,2]

	if ($isLeft) {
		### x-scale is scaled by 2 and y is scaled by 4/sqrt(3)
		@t1 = (3, 0); #(1.0*$scalex, 0*$scaley); ### [2,0]
		@t2 = (0, 2); #(-0.5*$scalex, sqrt(3.0)/2.0 *$scaley); ### [1,2]
	}

	my $ncrows = $ncmatrix->{"rows"};
	my $counter = 0;

	for (my $i = 0; $i < $lx; ++$i) {
		my $xv = $i*4;
		my $txv = 3*$i;
		#Sublattice A
		my $x0 = $xv;
		my $y0 = 0;
		my $tx0 = $txv;
		my $ty0 = 0;
		#Sublattice B
		my $sign = ($isLeft) ? -1 : 1;
		my $x1 = $xv + $sign;
		my $y1 = 1.0;
		my $tx1 = $txv + 0.5*$sign;
		my $ty1 = 1./$scaley;
		#Sublattice C
		my $x2 = $xv + 2*$sign;
		my $y2 = $y1;
		my $tx2 = $txv + 1.5*$sign;
		my $ty2 = $ty1;
		#Sublattice D
		my $x3 = $xv + 3*$sign;
		my $y3 = $y0;
		my $tx3 = $txv + 2*$sign;
		my $ty3 = $ty0;

		for (my $j = 0; $j < $ly; ++$j) {
			my $cxa = int($x0 + $j*$t2[0]);
			my $cxb = int($x1 + $j*$t2[0]);
			my $cxc = int($x2 + $j*$t2[0]);
			my $cxd = int($x3 + $j*$t2[0]);
			my $cya = int($y0 + $j*$t2[1]);
			my $cyb = int($y1 + $j*$t2[1]);
			my $cyc = int($y2 + $j*$t2[1]);
			my $cyd = int($y3 + $j*$t2[1]);

			my $txa = $tx0 + $j*$t2[0]/$scalex;
			my $txb = $tx1 + $j*$t2[0]/$scalex;
			my $txc = $tx2 + $j*$t2[0]/$scalex;
			my $txd = $tx3 + $j*$t2[0]/$scalex;
			my $tya = $ty0 + $j*$t2[1]/$scaley;
			my $tyb = $ty1 + $j*$t2[1]/$scaley;
			my $tyc = $ty2 + $j*$t2[1]/$scaley;
			my $tyd = $ty3 + $j*$t2[1]/$scaley;

			$tindx->[$counter] = $txa;
			$tindy->[$counter] = $tya;
			#print x0, y0, cxa, cya
			push @$ptsx, $cxa;
			push @$ptsy, $cya;
			$ncmatrix->{"data"}->[$cxa + $cya*$ncrows] = $counter++;

			$tindx->[$counter] = $txb;
			$tindy->[$counter] = $tyb;
			push @$ptsx, $cxb;
			push @$ptsy, $cyb;
			$ncmatrix->{"data"}->[$cxb + $cyb*$ncrows] = $counter++;

			$tindx->[$counter] = $txc;
			$tindy->[$counter] = $tyc;
			push @$ptsx, $cxc;
			push @$ptsy, $cyc;
			$ncmatrix->{"data"}->[$cxc + $cyc*$ncrows] = $counter++;
			
			$tindx->[$counter] = $txd;
			$tindy->[$counter] = $tyd;
			push @$ptsx, $cxd;
			push @$ptsy, $cyd;
			$ncmatrix->{"data"}->[$cxd + $cyd*$ncrows] = $counter++;
		}
	}
}



1;


