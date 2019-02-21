#!/usr/bin/perl

use strict;
use warnings;
use Math::Trig;

package OmegaUtils;

my $pi = Math::Trig::pi;

sub getLabels
{
	my ($hptr,$file) = @_;

	open(FILE,$file) or die "$0: Cannot open $file : $!\n";
	while (<FILE>) {
		chomp;
		foreach my $key (keys %$hptr) {
			if (/$key[= ]([^ ]+)/) {
				${$hptr->{$key}} = $1;
			}
		}
	}

	close(FILE);

	foreach my $key (keys %$hptr) {
		my $x = ${$hptr->{$key}};
		defined($x) or die "$0: Could not find $key in $file\n";
	}
}

sub printGnuplot
{
	my ($ptr, $geometry, $isPeriodic, $zeroAtCenter) = @_;

	my ($factor, $fileIndices, $leg) = getGeometryDetails($geometry);

	foreach my $fileIndex (@$fileIndices) {
		my $outFile = "outSpectrum$fileIndex.gnuplot";
		my $outFile2 = "outSpectrum$fileIndex.pgfplots";
		open(FOUT, ">", "$outFile") or die "$0: Cannot write to $outFile : $!\n";
		open(FOUT2, ">", "$outFile2") or die "$0: Cannot write to $outFile2 : $!\n";

		for my $omega (sort {$a <=> $b} keys %$ptr) {
			my $aptr = $ptr->{$omega};
			my $nks = scalar(@$aptr) - 1;
			my $numberOfQs = int($factor*$nks);
			my $centerShift = ($numberOfQs & 1) ? ($numberOfQs - 1)/2 : $numberOfQs/2;
			$centerShift = 0 unless ($zeroAtCenter);
			for (my $m2 = 0; $m2 < $numberOfQs; ++$m2) {
				my $m = $m2 - $centerShift;
				$m += $numberOfQs if ($m < 0);

				my $q = getQ($m2 - $centerShift, $numberOfQs, $isPeriodic);
				my $realPart = $aptr->[2*$m+1+2*$fileIndex*$numberOfQs];
				my $imagPart = $aptr->[2*$m+2+2*$fileIndex*$numberOfQs];
				print FOUT "$q $omega $realPart $imagPart\n";
				print FOUT2 "$q $omega $imagPart\n";
			}

			print FOUT2 "\n";
		}

		close(FOUT);
		close(FOUT2);
		print "$0: Written $outFile\n";
		print "$0: Written $outFile2\n";
	}
}

sub printOffsetPlots
{
	my ($ext, $ptr, $geometry, $isPeriodic, $zeroAtCenter) = @_;

	my ($factor, $fileIndices, $leg) = getGeometryDetails($geometry);

	foreach my $fileIndex (@$fileIndices) {
		my $offset = 0.7*findMaxVertical($ptr, $factor, $fileIndex);
		my $outFile = "outSpectrum$fileIndex.$ext";
		open(FOUT, ">", "$outFile") or die "$0: Cannot write to $outFile : $!\n";

		for my $omega (sort {$a <=> $b} keys %$ptr) {
			my $aptr = $ptr->{$omega};
			my $nks = scalar(@$aptr) - 1;
			my $numberOfQs = int($factor*$nks);
			my $centerShift = ($numberOfQs & 1) ? ($numberOfQs - 1)/2 : $numberOfQs/2;
			$centerShift = 0 unless ($zeroAtCenter);
			print FOUT "$omega ";
			for (my $m2 = 0; $m2 < $numberOfQs; ++$m2) {
				my $m = $m2 - $centerShift;
				$m += $numberOfQs if ($m < 0);

				my $q = getQ($m2 - $centerShift, $numberOfQs, $isPeriodic);
				#my $realPart = $aptr->[2*$m+1+2*$fileIndex*$numberOfQs];
				my $imagPart = $aptr->[2*$m+2+2*$fileIndex*$numberOfQs];
				my $val = $imagPart + $m2*$offset;
				print FOUT "$val ";
			}

			print FOUT "\n";
		}

		close(FOUT);
		print "$0: Written $outFile\n";
	}
}

sub getGeometryDetails
{
	my ($geometry, $my) = @_;
	my $factor = 0;
	my @fileIndices=(0);
	my $leg = 1;
	my $name = $geometry->{"name"};
	my $subname = $geometry->{"subname"};
	if ($name eq "chain") {
		$factor = 0.5;
		die "$0: Chain does not have ky != 0\n" if (defined($my) and $my != 0)
	} elsif ($name eq "ladder" || $subname eq "average") {
		$leg = $geometry->{"leg"};
		$factor = 0.25;
		@fileIndices=(0, 1) if ($leg == 2 || $subname eq "average");
	} else {
		die "$0: Unknown geometry $name\n";
	}

	return ($factor, \@fileIndices, $leg);
}



sub findMaxVertical
{
	my ($ptr, $factor, $fileIndex) = @_;

	my $max = 0;
	for my $omega (sort {$a <=> $b} keys %$ptr) {
		my $aptr = $ptr->{$omega};
		my $nks = scalar(@$aptr) - 1;
		my $numberOfQs = int($factor*$nks);
		for (my $m = 0; $m < $numberOfQs; ++$m) {
			my $imagPart = abs($aptr->[2*$m + 2 + 2*$fileIndex*$numberOfQs]);
			$max = $imagPart if ($max < $imagPart);
		}
	}

	return $max;
}

sub fourier
{
	my ($f, $v, $geometry, $hptr) = @_;

	my $subname = $geometry->{"subname"};

	if ($subname eq "average") {
		return fourierLadderAverage($f, $v, $geometry->{"leg"}, $hptr);
	}

	my $name = $geometry->{"name"};
	if ($name eq "chain") {
		return fourierChain($f, $v, $hptr);
	}

	if ($name eq "ladder") {
		return fourierLadder($f, $v, $geometry->{"leg"}, $hptr);
	}

	if ($name eq "LongRange" || $name eq "General") {

		my $orbitals = $hptr->{"Orbitals"};

		if ($geometry->{"subname"} eq "chain") {
			if ($orbitals == 1) {
				return fourierChain($f, $v, $hptr);
			} elsif ($orbitals == 2) {
                                return fourierChain2orb($f, $v, $hptr);
			}
	}

	if ($geometry->{"subname"} eq "ladder") {
                        if ($orbitals == 1) {
                                return fourierLadder($f, $v, $geometry->{"leg"}, $hptr);
                        } elsif ($orbitals == 2) {
                                return fourierLadder2orb($f, $v, $hptr);
                        }
		}
        }

	die "$0: ft: undefined geometry ".$geometry->{"name"}."\n";
}

sub fourierChain
{
	my ($f, $v, $hptr) = @_;
	my $n = scalar(@$v);
	my $mMax = $hptr->{"mMax"};
	my $isPeriodic = $hptr->{"isPeriodic"};
	my $centralSite = $hptr->{"centralSite"};
	my $numberOfQs = (defined($mMax)) ? $mMax : $n;
	for (my $m = 0; $m < $numberOfQs; ++$m) {
		my @sum = (0,0);
		my $q = getQ($m, $numberOfQs, $isPeriodic);
		for (my $i = 0; $i < $n; $i++) {
			my $ptr = $v->[$i];
			my @temp = @$ptr;
			my $arg = $q*($i-$centralSite);
			my $carg = cos($arg);
			my $sarg = sin($q*($i + 1));
			my $cOrSarg = ($isPeriodic) ? $carg : $sarg;
			$sum[0] += $temp[0]*$cOrSarg;
			$sum[1] += $temp[1]*$cOrSarg;
		}

		$f->[$m] = \@sum;
	}
}

sub fourierLadder
{
	my ($f, $v, $leg, $hptr) = @_;
	my $n = scalar(@$v);
	my $mMax = $hptr->{"mMax"};
	my $isPeriodic = $hptr->{"isPeriodic"};
	my $centralSite = $hptr->{"centralSite"};
	my $numberOfQs = (defined($mMax)) ? $mMax : int($n/$leg);
	for (my $m = 0; $m < $numberOfQs; ++$m) {
		my $q = getQ($m, $numberOfQs, $isPeriodic);

		my @fPerLeg;
		my @sumKy0 = (0, 0);
		for (my $ll = 0; $ll < $leg; ++$ll) {
			my @f = fourierF($v, $q, $ll, $leg, $centralSite);
			$fPerLeg[$ll] = \@f;
			$sumKy0[0] += $f[0];
			$sumKy0[1] += $f[1];
		}

		$f->[$m] = \@sumKy0;

		next if ($leg != 2);

		for (my $x = 1; $x < 2; ++$x) {
			my $sign = 1-2*$x;
			my $realPart = $fPerLeg[0]->[0] + $sign*$fPerLeg[1]->[0];
			my $imagPart = $fPerLeg[0]->[1] + $sign*$fPerLeg[1]->[1];
			my @sum = ($realPart,$imagPart);
			$f->[$m+$numberOfQs*$x] = \@sum;
		}
	}
}

sub fourierF
{
	my ($v, $q, $ll, $leg, $centralSite) = @_;
	my $n = scalar(@$v);
	my @sum;
	for (my $i = $ll; $i < $n; $i += $leg) {
		my $ptr = $v->[$i];
		my @temp = @$ptr;
		my $arg = $q*distanceLadder($i, $centralSite, $ll, $leg);
		# FIXME THINK ABOUT OPEN BC
		my $carg = cos($arg);
		$sum[0] += $temp[0]*$carg;
		$sum[1] += $temp[1]*$carg;
	}

	return @sum;
}

sub distanceLadder
{
	my ($ind, $jnd, $ll, $leg) = @_;
	return ($ind - $ll - $jnd)/$leg;
}

sub fourierChain2orb
{
        my ($f, $v, $hptr) = @_;
	my $mMax = $hptr->{"mMax"};
	my $isPeriodic = $hptr->{"isPeriodic"};
        my $n = int(0.25*scalar(@$v));
        my $numberOfQs = (defined($mMax)) ? $mMax : $n;
	my $cSite = $n/2-1;
        for (my $m = 0; $m < $numberOfQs; ++$m) {
                my @sum = (0,0);
                my $q = getQ($m, $numberOfQs, $isPeriodic);

# orb A
		for (my $i = 0; $i < $n; $i++) {
			my @temp = ($v->[4*$i],$v->[4*$i+1]);
			my $arg = $q*($i-$cSite);
			my $carg = cos($arg); 
			$sum[0] += $temp[0]*$carg;
			$sum[1] += $temp[1]*$carg;
		}
# orb B
		for (my $i = 0; $i < $n; $i++) {
			my @temp = ($v->[4*$i+2],$v->[4*$i+3]);
			my $arg = $q*($i-$cSite);
			my $carg = cos($arg);
			$sum[0] += $temp[0]*$carg;
			$sum[1] += $temp[1]*$carg;
		}

		$f->[$m] = \@sum;
	}
}

sub fourierLadder2orb
{
        my ($f, $v, $leg, $hptr) = @_;
	my $mMax = $hptr->{"mMax"};
	my $isPeriodic = $hptr->{"isPeriodic"};
	my $centralSite = $hptr->{"centralSite"};
        my $n = int(0.125*scalar(@$v));
        my $numberOfQs = (defined($mMax)) ? $mMax : $n;
        for (my $m = 0; $m < $numberOfQs; ++$m) {
                my $q = getQ($m,$numberOfQs,$isPeriodic);
                my @f0 = fourierF0test($v, $q, $centralSite);
                my @f1 = fourierF1test($v, $q, $centralSite);
                for (my $x = 0; $x < 2; ++$x) {
                        my $sign = 1-2*$x;
                        my $realPart = $f0[0] + $sign*$f1[0];
                        my $imagPart = $f0[1] + $sign*$f1[1];
                        my @sum = ($realPart,$imagPart);
                        $f->[$m+$numberOfQs*$x] = \@sum;
                }
        }
}


sub fourierF0test
{
        my ($v,$q) = @_;
        my $n = int(0.25*scalar(@$v));
	my $cSite = $n/2-2;
        my @sum;
# orb A
        for (my $i = 0; $i < $n; $i+=2) {
                my @temp = ($v->[4*$i],$v->[4*$i+1]);
                my $arg = $q*($i-$cSite)*0.5;
                my $carg = cos($arg);
                $sum[0] += $temp[0]*$carg;
                $sum[1] += $temp[1]*$carg;
        }

# orb B
        for (my $i = 0; $i < $n; $i+=2) {
                my @temp = ($v->[4*$i+2],$v->[4*$i+3]);
                my $arg = $q*($i-$cSite)*0.5;
                my $carg = cos($arg);
                $sum[0] += $temp[0]*$carg;
                $sum[1] += $temp[1]*$carg;
        }

        return @sum;
}


sub fourierF1
{
	my ($v, $q, $centralSite) = @_;
	my $n = int(0.5*scalar(@$v));
	my @sum;
	for (my $i = 1; $i < $n; $i+=2) {
		my @temp = ($v->[2*$i],$v->[2*$i+1]);
		my $arg = $q*distanceLadder($i,$centralSite,0,2);
		my $carg = cos($arg);
		$sum[0] += $temp[0]*$carg;
		$sum[1] += $temp[1]*$carg;
	}

	return @sum;
}

sub fourierF1test
{
        my ($v,$q) = @_;
        my $n = int(0.25*scalar(@$v));
        my @sum;
        my $cSite = $n/2-2;
# orb A
        for (my $i = 1; $i < $n; $i+=2) {
                my @temp = ($v->[4*$i],$v->[4*$i+1]);
                my $arg = $q*distanceLadder($i,$cSite,0,2);
                my $carg = cos($arg);
                $sum[0] += $temp[0]*$carg;
                $sum[1] += $temp[1]*$carg;
        }
# orb B
        for (my $i = 1; $i < $n; $i+=2) {
                my @temp = ($v->[4*$i+2],$v->[4*$i+3]);
                my $arg = $q*distanceLadder($i,$cSite,0,2);
                my $carg = cos($arg);
                $sum[0] += $temp[0]*$carg;
                $sum[1] += $temp[1]*$carg;
        }

        return @sum;
}

sub fourierLadderAverage
{
	my ($f, $v, $leg, $hptr) = @_;
	my $n = scalar(@$v);
	my $lx = int($n/$leg);
	my $legSmall = 2;
	my $total = int($leg/$legSmall);
	my $centralSite = $hptr->{"centralSite"};
	my $ll = ($centralSite == int($n/2)) ? 0 : 1;
	die "$0: Wrong central site $centralSite\n" if ($centralSite != int($n/2) + $ll*2);

	# prepare partialV

	my @partialV = fourierLadderAverageInput($v, $ll, $leg, $hptr);
	print STDERR "$0: partialV for ll=$ll with ".scalar(@partialV)." entries\n";
	my %modifiedHptr = %$hptr;
	$modifiedHptr{"centralSite"} = $lx;
	my @partialF;
	fourierLadder($f, \@partialV, 2, \%modifiedHptr);
}

sub fourierLadderAverageInput
{
	my ($v, $ll, $leg, $hptr) = @_;
	my $n = scalar(@$v);
	my $lx = int($n/$leg);
	my $legSmall = 2;
	my @partialV;
	for (my $x = 0; $x < $lx; ++$x) {
		for (my $y = 0; $y < $legSmall; ++$y) {
			my $oldY = $y + 2*$ll;
			$partialV[$y + $x*$legSmall] = $v->[$oldY + $x*$leg];
		}
	}

	return @partialV;
}

sub getQ
{
	my ($m, $n, $isPeriodic) = @_;
	return ($isPeriodic) ? 2.0*$pi*$m/$n : ($m + 1)*$pi/($n+1.0);
}

1;

