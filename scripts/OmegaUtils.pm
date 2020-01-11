#!/usr/bin/perl

use strict;
use warnings;
use Math::Trig;
use Honeycomb;

package OmegaUtils;

my $pi = Math::Trig::pi;

sub isAinur
{
	my ($file) = @_;
	open(FILE, "<", "$file") or die "$0: Cannot open $file : $!\n";
	$_ = <FILE>;
	close(FILE);
	chomp;
	return $_ eq "##Ainur1.0";
}

sub getLabels
{
	my ($hptr,$file) = @_;

	open(FILE,$file) or die "$0: Cannot open $file : $!\n";
	while (<FILE>) {
		chomp;
		foreach my $key (keys %$hptr) {
			if (/$key[= ]([^ ]+)/) {
				my $newVal = $1;
				my $prev = ${$hptr->{$key}};
				if ($prev && $prev ne $newVal) {
					print STDERR "Already a previous value for $key of $prev\n";
					print "New value is $newVal\n";
					print "To take new value press ENTER. Or enter value ";
					$_ = <STDIN>;
					chomp;
					if ($_) {
						${$hptr->{$key}} = $_;
						print STDERR "$0: Value for $key is ".${$hptr->{$key}}."\n";
						next;
					}
				}

				${$hptr->{$key}} = $newVal;
				print STDERR "$0: Value for $key is ".${$hptr->{$key}}."\n";
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
	my ($inFile, $geometry, $isPeriodic, $zeroAtCenter, $nonNegativeOnly) = @_;
	my $hasPrinted = 0;
	open(FIN, "<", "$inFile") or die "$0: Cannot open $inFile : $!\n";

	my %h;
	my $npoints;
	while (<FIN>) {
		my @temp = split;
		my $n = scalar(@temp);
		if ($n < 1) {
			print STDERR "$0: line $. in $inFile is empty, skipping\n";
			next;
		}

		print STDERR "$0: Columns $n in $inFile\n" if (!$hasPrinted);
		$hasPrinted = 1;

		my $omega = $temp[0];
		my $tmpPoints = scalar(@temp);
		defined($npoints) or $npoints = $tmpPoints;
		if ($npoints != $tmpPoints) {
			print "$0: Discared values for omega=$omega; found $tmpPoints, expected $npoints\n";
			next;
		}

		$h{$omega} = \@temp;
	}

	close(FIN);

	printGnuplotFromHash(\%h, $geometry, $isPeriodic, $zeroAtCenter, $nonNegativeOnly);
}

sub printGnuplotFromHash
{
	my ($ptr, $geometry, $isPeriodic, $zeroAtCenter, $nonNegativeOnly) = @_;

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
			my @array;
			if ($geometry->{"subname"} =~ /^Honeycomb/) {
				my $type = $geometry->{"subname"};
				@array = Honeycomb::fillQvalues($ptr, $type);
			} else {
				@array = fillQvalues($numberOfQs, $isPeriodic, $zeroAtCenter);
			}

			foreach my $ptr2 (@array) {
				my ($m, $q) = ($ptr2->{"m"}, $ptr2->{"q"});

				my $realPart = $aptr->[2*$m+1+2*$fileIndex*$numberOfQs];
				my $imagPart = $aptr->[2*$m+2+2*$fileIndex*$numberOfQs];
				$imagPart = 0 if ($nonNegativeOnly and $imagPart < 0);
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

sub fillQvalues
{
	my ($numberOfQs, $isPeriodic, $zeroAtCenter) = @_;
	my $centerShift = ($numberOfQs & 1) ? ($numberOfQs - 1)/2 : $numberOfQs/2;
	$centerShift = 0 unless ($zeroAtCenter);
	my @array;
	for (my $m2 = 0; $m2 < $numberOfQs; ++$m2) {
		my $m = $m2 - $centerShift;
		$m += $numberOfQs if ($m < 0);
		my $q = getQ($m2 - $centerShift, $numberOfQs, $isPeriodic);
		my $temp = { "m" => $m, "q" => $q};
		push @array, $temp;
	}

	return @array;
}

sub printOffsetPlots
{
	my ($ext, $ptr, $geometry, $isPeriodic, $zeroAtCenter) = @_;

	if ($geometry->{"subname"} =~ /^Honeycomb/) {
		print STDERR "$0: WARNING: Honeycomb: printOffsetPlots not supported yet\n";
		return;
	}

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
	} elsif ($subname eq "GrandCanonical" and $name eq "ladder") {
		$factor = 0.5;
		die "$0: Chain does not have ky != 0\n" if (defined($my) and $my != 0);
	} elsif ($name eq "ladder" || $subname eq "average") {
		$leg = $geometry->{"leg"};
		$factor = 0.25;
		@fileIndices=(0, 1) if ($leg == 2 || $subname eq "average");
	} elsif ($geometry->{"subname"} =~ /^Honeycomb/) {
		return (1, \@fileIndices, $leg);
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

	if ($name eq "ladder" and $subname eq "GrandCanonical") {
		return fourierChainGC($f, $v, $hptr);
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

		if ($geometry->{"subname"} =~ /^Honeycomb/) {
			my $type = $geometry->{"subname"};
			$type =~ s/^Honeycomb//;
			return fourierHoneycomb($f, $v, $hptr, $type);
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
	my @centralSites = ($centralSite);

	if (!$isPeriodic) {
		my $b = ($centralSite == int($n/2));
		if (!$b && ($centralSite != int($n/2) - 1)) {
			die "$0: Chain of $n sites, but central site is $centralSite, makes no sense!?\n";
		}

		# FIXME: DOES NOT WORK, CHECK FORMULA BELOW
		if ($hptr->{"multicenter"}) {
			my $otherCenter = ($b) ? $centralSite - 1 : $centralSite + 1;
			push @centralSites, $otherCenter;
		}
	}

	my $numberOfQs = (defined($mMax)) ? $mMax : $n;
	for (my $m = 0; $m < $numberOfQs; ++$m) {
		my @sum = (0,0);
		my $q = getQ($m, $numberOfQs, $isPeriodic);
		foreach my $cSite (@centralSites) {
			for (my $i = 0; $i < $n; $i++) {
				my $ptr = $v->[$i];
				my @temp = @$ptr;
				my $arg = $q*($i-$cSite);
				my $carg = cos($arg);
				my $sarg = sin($q*($i + 1))*sin($q*($cSite + 1));
				my $cOrSarg = ($isPeriodic) ? $carg : $sarg;
				$sum[0] += $temp[0]*$cOrSarg;
				$sum[1] += $temp[1]*$cOrSarg;
			}
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
			my @f = fourierF($v, $q, $ll, $leg, $centralSite, $isPeriodic);
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
	my ($v, $q, $ll, $leg, $centralSite, $isPeriodic) = @_;
	my $n = scalar(@$v);
	my @sum;
	for (my $i = $ll; $i < $n; $i += $leg) {
		my $ptr = $v->[$i];
		my @temp = @$ptr;
		my $arg = $q*distanceLadder($i, $centralSite, $ll, $leg);
		my $carg = cos($arg);
		my $sarg = sin($q*($i + 1))*sin($q*($centralSite + 1));
		my $cOrSarg = ($isPeriodic) ? $carg : $sarg;
		$sum[0] += $temp[0]*$cOrSarg;
		$sum[1] += $temp[1]*$cOrSarg;
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

		for (my $orb = 0; $orb < 2; ++$orb) {
			for (my $i = 0; $i < $n; $i++) {
				my @temp = ($v->[4*$i + 2*$orb],$v->[4*$i + 2*$orb + 1]);
				my $arg = $q*($i-$cSite);
				my $carg = cos($arg);
				my $sarg = sin($q*($i + 1))*sin($q*($cSite + 1));
				my $cOrSarg = ($isPeriodic) ? $carg : $sarg;
				$sum[0] += $temp[0]*$cOrSarg;
				$sum[1] += $temp[1]*$cOrSarg;
			}
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
		my @f0 = fourierF0test($v, $q, $centralSite, $isPeriodic);
		my @f1 = fourierF1test($v, $q, $centralSite, $isPeriodic);
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
	my ($v, $q, $cSite, $isPeriodic) = @_;
	my $n = int(0.25*scalar(@$v));
	my @sum;
	for (my $orb = 0; $orb < 2; ++$orb) {
		for (my $i = 0; $i < $n; $i += 2) {
			my @temp = ($v->[4*$i + 2*$orb],$v->[4*$i + 2*$orb + 1]);
			my $arg = $q*($i-$cSite)*0.5;
			my $carg = cos($arg);
			my $sarg = sin($q*($i + 1))*sin($q*($cSite + 1));
			my $cOrSarg = ($isPeriodic) ? $carg : $sarg;
			$sum[0] += $temp[0]*$cOrSarg;
			$sum[1] += $temp[1]*$cOrSarg;
		}
	}

	return @sum;
}

sub fourierF1
{
	my ($v, $q, $cSite, $isPeriodic) = @_;
	my $n = int(0.5*scalar(@$v));
	my @sum;
	for (my $i = 1; $i < $n; $i+=2) {
		my @temp = ($v->[2*$i],$v->[2*$i+1]);
		my $arg = $q*distanceLadder($i,$cSite,0,2);
		my $carg = cos($arg);
		my $sarg = sin($q*($i + 1))*sin($q*($cSite + 1));
		my $cOrSarg = ($isPeriodic) ? $carg : $sarg;
		$sum[0] += $temp[0]*$cOrSarg;
		$sum[1] += $temp[1]*$cOrSarg;
	}

	return @sum;
}

sub fourierF1test
{
	my ($v, $q, $cSite, $isPeriodic) = @_;
	my $n = int(0.25*scalar(@$v));
	my @sum;
	for (my $orb = 0; $orb < 2; ++$orb) {
		for (my $i = 1; $i < $n; $i += 2) {
			my @temp = ($v->[4*$i + 2*$orb],$v->[4*$i + 2*$orb + 1]);
			my $arg = $q*distanceLadder($i,$cSite,0,2);
			my $carg = cos($arg);
			my $sarg = sin($q*($i + 1))*sin($q*($cSite + 1));
			my $cOrSarg = ($isPeriodic) ? $carg : $sarg;
			$sum[0] += $temp[0]*$cOrSarg;
			$sum[1] += $temp[1]*$cOrSarg;
		}
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

sub fourierHoneycomb
{
	my ($f, $v, $hptr, $type) = @_;
	my $lx = $hptr->{"#Lx"};
	my $ly = $hptr->{"#Ly"};
	my ($M1, $M2) = (2*$lx, $ly);
	my $n = scalar(@$v);
	my (@tindx, @tindy);
	Honeycomb::honeySpace(\@tindx, \@tindy, $n, $hptr, $type);
	for (my $m1 = 0; $m1 < $M1; ++$m1) { # loop over momenta
		for (my $m2 = 0; $m2 < $M2; ++$m2) { # loop over momenta

			# valid ($qx, $qy)
			my ($qx, $qy) = Honeycomb::honeyGetQ($m1, $m2, $lx, $ly, $type);
			my @sum = (0, 0);
			for (my $i = 0; $i < $n; ++$i) { # loop over space
				my $ptr = $v->[$i];
				my @temp = @$ptr;
				my @fourierFactor = honeyFourierFactor($i, $qx, $qy, \@tindx, \@tindy, $hptr);
				$sum[0] += $fourierFactor[0]*$temp[1] + $fourierFactor[1]*$temp[0]; # imaginary part
				$sum[1] += $fourierFactor[0]*$temp[0] - $fourierFactor[1]*$temp[1]; # real part
			}

			$f->[$m1 + $m2*$M1] = \@sum;
		}
	}
}

sub honeyFourierFactor
{
	my ($i, $qx, $qy, $tindx, $tindy, $hptr) = @_;
	# get (rx, ry) and (cx, cy)
	my $indexForCenter = $hptr->{"centralSite"};
	my ($rx, $ry) = ($tindx->[$i], $tindy->[$i]);
	my ($cx, $cy) = ($tindx->[$indexForCenter], $tindy->[$indexForCenter]);
	my $arg = $qx*($rx - $cx) + $qy*($ry - $cy);
	return (cos($arg), sin($arg));
}

sub fourierChainGC
{
	my ($f, $v, $hptr) = @_;
	my $n = scalar(@$v);

	my $mMax = $hptr->{"mMax"};
	my $isPeriodic = $hptr->{"isPeriodic"};
	my $centralSite = $hptr->{"centralSite"};
	my $nOver2 = int($n/2);

	die "$0: FATAL: ChainGC central site is odd\n" if ($centralSite & 1);

	if (!$isPeriodic) {
		my $b = ($centralSite != $nOver2);
		if ($b && $centralSite != $nOver2 - 2) {
			die "$0: FATAL ChainGC: wrong central site $centralSite\n";
		}
	}
	
	my $cSite = int($centralSite/2);
	my $numberOfQs = (defined($mMax)) ? $mMax : $nOver2;
	for (my $m = 0; $m < $numberOfQs; ++$m) {
		my @sum = (0,0);
		my $q = getQ($m, $numberOfQs, $isPeriodic);
		for (my $ii = 0; $ii < $n; $ii += 2) {
			my $i = int($ii/2);
			my $ptr = $v->[$ii];
			my @temp = @$ptr;
			my $arg = $q*($i - $cSite);
			my $carg = cos($arg);
			my $sarg = sin($q*($i + 1))*sin($q*($cSite + 1));
			my $cOrSarg = ($isPeriodic) ? $carg : $sarg;
			$sum[0] += $temp[0]*$cOrSarg;
			$sum[1] += $temp[1]*$cOrSarg;
		}

		$f->[$m] = \@sum;
	}
}

sub writeFourier
{
	my ($array, $f, $geometry) = @_;
	my $subname = $geometry->{"subname"};
	my $isPeriodic = $geometry->{"isPeriodic"};

	if ($geometry->{"name"} eq "chain" || $subname eq "GrandCanonical") {
		return writeFourierChain($array,$f, $isPeriodic);
	}

	if ($geometry->{"name"} eq "ladder" || $subname eq "average") {
		return writeFourierLadder($array, $f);
	}

	die "$0: writeFourier: undefined geometry ".$geometry->{"name"}."\n";
}

sub writeFourierChain
{
	my ($array, $f, $isPeriodic) = @_;

	my $n = scalar(@$f);
	for (my $m = 0; $m < $n; ++$m) {
		my $q = getQ($m, $n, $isPeriodic);
		my $ptr = $f->[$m];
		my @temp = @$ptr;
		$array->[$m] = [$q, $temp[0], $temp[1]];
	}
}

sub writeFourierLadder
{
	my ($array, $f) = @_;

	my $n = scalar(@$f);
	for (my $m = 0; $m < $n; ++$m) {
		my $ptr = $f->[$m];
		my @temp = @$ptr;
		my @temp2 = ($m);
		push @temp2, @temp;
		$array->[$m] = \@temp2;
	}
}

sub getQ
{
	my ($m, $n, $isPeriodic) = @_;
	return ($isPeriodic) ? 2.0*$pi*$m/$n : ($m + 1)*$pi/($n+1.0);
}

1;

