#!/usr/bin/perl
# All credit goes to Nirav P.
# All errors go to G.A.

use strict;
use warnings;
use utf8;
use List::Util qw( min max );

my ($templateInput, $templateTex) = @ARGV;
defined($templateInput) or die "USAGE: $0 templateInput [templateTex]\n";

my @letters = ("x", "y", "z");
my @magnetic = ("BxBxBx", "ByByBy", "BzBzBz");
my @Kit = ("KxKxKx", "KyKyKy", "KzKzKz");
my @Gammas = ("GammaX", "GammaY", "GammaZ");
my @GammaPrimes = ("GammaPrimeX", "GammaPrimeY", "GammaPrimeZ");
my $options = "periodicy";
my ($bx, $by, $bz) = @magnetic;
my ($kxkx, $kyky, $kzkz) = @Kit;
my ($gammaX, $gammaY, $gammaZ) = @Gammas;
my ($gammaPrimeX, $gammaPrimeY, $gammaPrimeZ) = @GammaPrimes;
my $printcut = 20;
my ($lx, $ly);
my ($jx, $jy, $jz) = (0,0,0);
my ($j3x, $j3y, $j3z) = (0,0,0);
my $model = "UNKNOWN";
my $dashed = 1;
my $withNumbers = 1;
my $hopt = {"Model" => \$model,
            "#lx" => \$lx,
            "#ly" => \$ly,
            "#options" => \$options,
            "#bx" => \$bx,
            "#by" => \$by,
            "#bz" => \$bz,
	    "#kxkx" => \$kxkx,
            "#kyky" => \$kyky,
            "#kzkz" => \$kzkz,
            "#jx" => \$jx,
	    "#jy" => \$jy,
            "#jz" => \$jz,
            "#j3x" => \$j3x,
	    "#j3y" => \$j3y,
            "#j3z" => \$j3z,
            "#gammax" => \$gammaX,
            "#gammay" => \$gammaY,
            "#gammaz" => \$gammaZ,
            "#gammaPrimeX" => \$gammaPrimeX,
            "#gammaPrimeY" => \$gammaPrimeY,
            "#gammaPrimeZ" => \$gammaPrimeZ,
            "#printcut" => \$printcut,
            "#dashed" => \$dashed,
            "#withNumbers" => \$withNumbers};

getLabels($hopt, $templateInput);
 
@magnetic = ($bx, $by, $bz);
@Kit = ($kxkx, $kyky, $kzkz);
@Gammas = ($gammaX, $gammaY, $gammaZ);
@GammaPrimes = ($gammaPrimeX, $gammaPrimeY, $gammaPrimeZ);
my @heisenbergJ = ($jx, $jy, $jz);
my @heisenbergJ3 = ($j3x, $j3y, $j3z);

my $isHash = procOptions($options);
my $isArmchairX = $isHash->{"isArmchairX"};
my $n = $lx*$ly*2;
if ($isArmchairX){
	$n = $lx*$ly*4;
}
my $replaceMagnetic;
for (my $dir = 0; $dir < 3; ++$dir) {
	$replaceMagnetic .= getVector("MagneticField" . uc($letters[$dir]), $n, $magnetic[$dir]);
	$replaceMagnetic .= "\n";
}

my ($replaceConnections, $plot) = getConnectionsAndPlot($lx, $ly, $options);

my $nOfTerms = 3;
if ($model eq "Kitaev") {
	;
} elsif ($model eq "KitaevExtended") {
	$nOfTerms = 5;
} elsif ($model eq "KitaevWithGammas") {
	$nOfTerms = 9;
} else {
	die "$0: Model $model not supported\n";
}

my $params = {"n" => $n, 
              "replaceMagnetic" => $replaceMagnetic,
              "replaceConnections" => $replaceConnections,
              "halfNminusOne" => ($n/2 - 1),
              "nMinusTwo" => ($n - 2),
              "nOfTerms" => $nOfTerms};
createInput("test.inp", $params, $templateInput);

my $params2 = {"plot" => $plot};

if (defined($templateTex)) {
	createInput("test.tex", $params2, $templateTex);
} else {
	print STDERR "$0: No templateTex given, no tex created\n";
}

sub createInput
{
	my ($file, $params, $templateInput) = @_;

	open(FOUT, ">", "$file") or die "$0: Cannot write to $file\n";

	open(FILE, "<", "$templateInput") or die "$0: Cannot open $templateInput: $!\n";

	while(<FILE>) {
		next if (/^#/);
		if (/\$([a-zA-Z0-9\[\]]+)/) {
				my $name = $1;
				my $str = "\$"."params->{\"$name\"}";
				my $val = eval "$str";
				defined($val) or die "$0: Undefined substitution for $name\n";
				s/\$\Q$name/$val/g;
		}
		print FOUT;
	}

	close(FILE);
	close(FOUT);
	print STDERR "$0: File $file has been written\n";

	return $file;
}

sub getConnectionsAndPlot
{
	my ($lx, $ly, $options) = @_;
	my $n = 2*$lx*$ly;
	my $gconnections = "";
	my $connections = "";
	my @n1neigh;
	my @n3neigh;
	my @indx;
	my @indy;

	my $isHash = procOptions($options);
	my $isArmchairX = $isHash->{"isArmchairX"};
	if ($isArmchairX){
		$n = $lx*$ly*4;
	}

	my $plot = honeycomb(\@n1neigh, \@n3neigh, \@indx, \@indy, $lx, $ly, $isHash);

	for (my $dir = 0; $dir < 3; ++$dir) {
		my @kthisdir;
		my $total = $n*$n;
		for (my $i = 0; $i < $total; ++$i) {
			$kthisdir[$i] = 0;
		}


		for (my $i = 0; $i < $n; ++$i) {
			#Kitaev
			my $j = $n1neigh[$i + $dir*$n];
			if (defined($j) and $i < $j) {
				#my $jx = int($indx[$j]);
				#my $jy = int($indy[$j]);
				addSymmetric(\@kthisdir, $i, $j, $n, $Kit[$dir]);
			}

			#Heisenberg
			for (my $dir2 = 0; $dir2 < 3; ++$dir2) {
				my $jj = $n1neigh[$i + $dir2*$n];
				next unless (defined($jj) and $i < $jj);
				addSymmetric(\@kthisdir, $i, $jj, $n, $heisenbergJ[$dir2]);
			}

			#Third nearest Heisenberg
			for (my $dir2 = 0; $dir2 < 3; ++$dir2) {
				my $jj = $n3neigh[$i + $dir2*$n];
				next unless (defined($jj) and $i < $jj);
				addSymmetric(\@kthisdir, $i, $jj, $n, $heisenbergJ3[$dir2]);
			}
		}

		my $kthismatrix = {"data" => \@kthisdir, "rows" => $n, "cols" => $n};

		$connections .= getDmrgppMatrix($kthismatrix);
		$connections .= "\n" unless ($dir == 2);
	}

	for (my $dir = 0; $dir < 6; ++$dir) {
		my @gammathisdir;
		my $total = $n*$n;
		for (my $i = 0; $i < $total; ++$i) {
			$gammathisdir[$i] = 0;
		}

		my $dirGamma = int($dir/2);
		for (my $i = 0; $i < $n; ++$i) {
			#Gammas
			my $j = $n1neigh[$i + $dirGamma*$n];
			if (defined($j) and $i < $j) {
				addSymmetric(\@gammathisdir, $i, $j, $n, $Gammas[$dirGamma]);
			}

			#Gamma primes // connects S^dir0 * S^dir1 
			my ($dir0, $dir1) = findGammaDirections($dir);
			my $jj = $n1neigh[$i + $dir0*$n];
			if (defined($jj) and $i < $jj) {
				addSymmetric(\@gammathisdir, $i, $jj, $n, $GammaPrimes[$dir0]);
			}

			$jj = $n1neigh[$i + $dir1*$n];
			if (defined($jj) and $i < $jj) {
				addSymmetric(\@gammathisdir, $i, $jj, $n, $GammaPrimes[$dir1]);
			}
		}

		my $gthismatrix = {"data" => \@gammathisdir, "rows" => $n, "cols" => $n};
		my $someString = getDmrgppMatrix($gthismatrix);
		$gconnections .= "$someString";
		$gconnections .= "\n" unless ($dir == 5);
	}

	if ($model eq "KitaevWithGammas") {
		$connections .= $gconnections;
	}

	return ($connections, $plot);
}

sub findGammaDirections
{
	my ($dir) = @_;
	my $dirGamma = int($dir/2);
	my @dirs;
	for (my $d = 0; $d < 3; ++$d) {
		next if ($d == $dirGamma);
		push @dirs, $d;
	}

	my $n = scalar(@dirs);
	($n == 2) or die "$0: findGammaDirections failed for dir=$dir\n";
	if ($dir & 1) {
		my $tmp = $dirs[0];
		$dirs[0] = $dirs[1];
		$dirs[1] = $tmp;
	}

	return @dirs;
}

sub addSymmetric
{
	my ($a, $i, $j, $n, $value) = @_;
	my $tmp = $a->[$i + $j*$n];

	if ($tmp =~ /^[\d\.\-\+e]+$/ and $value =~ /^[\d\.\-\+e]+$/) {
		$a->[$i + $j*$n] += $value;
		$a->[$j + $i*$n] += $value;
	} else {
		my $extra = "";
		if ($tmp eq "0") {
			$a->[$i + $j*$n] = $a->[$j + $i*$n] = "";
		} else {
			$extra = " + ";
		}

		$a->[$i + $j*$n] .= "$extra$value";
		$a->[$j + $i*$n] .= "$extra$value";
	}
}

sub getDmrgppMatrix
{
	my ($m) = @_;
my $str= <<EOF;
DegreesOfFreedom=1
GeometryKind=LongRange
GeometryOptions=none
EOF
	return $str . getMatrix("Connectors", $m);
}

sub getVector
{
	my ($label, $total, $value) = @_;
	my $str = "$label $total\n";
	for (my $i = 0; $i < $total; ++$i) {
		$str .= "$value";
		if (($i + 1) % $printcut == 0) {
			$str .= "\n";
		} else {
			$str .= " " if ($i + 1 != $total);
		}
	}

	return "$str\n";
}

sub getMatrix
{
	my ($label, $m) = @_;
	my $rows = $m->{"rows"};
	my $cols = $m->{"cols"};
	my $total = $rows*$cols;
	my $str = "$label $rows $cols\n";
	for (my $i = 0; $i < $total; ++$i) {
		my $value = $m->{"data"}->[$i];
		defined($value) or die "$i\n";
		$str .= "$value ";
		if (($i + 1) % $cols == 0) {
			$str .= "\n";
		} else {
			$str .= " " if ($i + 1 != $total);
		}
	}

	return "$str\n";
}

sub procOptions
{
	my ($options) = @_;
	my $isPeriodicX = ($options =~ /periodicx/i);
	my $isPeriodicY = ($options =~ /periodicy/i);
	my $isArmchairX = ($options =~ /armchairx/i);
	my $isLeft = ($options =~ /left/i);
	return {"isPeriodicX" => $isPeriodicX,
	        "isPeriodicY" => $isPeriodicY,
	        "isArmchairX" => $isArmchairX,
	        "isLeft" => $isLeft};
}

sub honeycomb
{
	my ($n1neigh, $n3neigh, $indx, $indy, $lx, $ly, $isHash) = @_;

	my $isPeriodicX = $isHash->{"isPeriodicX"};
	my $isPeriodicY = $isHash->{"isPeriodicY"};
	my $isArmchairX = $isHash->{"isArmchairX"};
	my $isLeft = $isHash->{"isLeft"};
		
	my $n = $lx*$ly*2;
	my $ncrows = $lx*2 + $ly;
	my $nccols = $ly*2;
	if ($isArmchairX){
		$n = $lx*$ly*4;
		$ncrows = $lx*4  + $ly;
		$nccols = $ly*2;
	}
	my @nc;
	my @tindx;
	my @tindy;
	my @ptsx;
	my @ptsy;
	my $ncmatrix = {"data" => \@nc, "rows" => $ncrows, "cols" => $nccols};
	print "\nncrows=", $ncrows, ", nccols=", $nccols,"\n";
	if ($isArmchairX){
		honeycombArmchairSetSiteCoordinates(\@tindx, \@tindy, \@ptsx, \@ptsy, $ncmatrix, $lx, $ly, $isHash);
	} else {
		honeycombZigzagSetSiteCoordinates(\@tindx, \@tindy, \@ptsx, \@ptsy, $ncmatrix, $lx, $ly, $isHash);
	}


	#### ===============================================
	my $xmax = max @ptsx;
	my $ymax = max @ptsy;

	for (my $i = 0; $i < $n; ++$i) {
		$indx->[$i] = int($ptsx[$i]);
		$indy->[$i] = int($ptsy[$i]);
	}

	my $n1neighrows = $n;
	my $n1neighcols = 3;
	my $n3neighrows = $n;
	my $n3neighcols = 3;

	my $n1neighmatrix = {"data" => $n1neigh, "rows" => $n1neighrows, "cols" => $n1neighcols};
	my $n3neighmatrix = {"data" => $n3neigh, "rows" => $n3neighrows, "cols" => $n3neighcols};
	
	if ($isArmchairX){
		honeycombArmchairSetNeigh($n1neighmatrix, $lx, $ly, $isHash, $indx, $indy, $ncmatrix, $xmax, $ymax);
		honeycombArmchairSetThirdNeigh($n3neighmatrix, $lx, $ly, $isHash, $indx, $indy, $ncmatrix, $xmax, $ymax);
		return honeycombArmchairPlot($indx, $indy, \@tindx, \@tindy, $n1neighmatrix, $ncmatrix, $isHash)
	} else {
		honeycombZigzagSetNeigh($n1neighmatrix, $lx, $ly, $isHash, $indx, $indy, $ncmatrix, $xmax, $ymax);
		honeycombZigzagSetThirdNeigh($n3neighmatrix, $lx, $ly, $isHash, $indx, $indy, $ncmatrix, $xmax, $ymax);	
		return honeycombZigzagPlot($indx, $indy, \@tindx, \@tindy, $n1neighmatrix, $ncmatrix, $isHash)
	}
	
}

sub honeycombZigzagSetNeigh
{
	my ($n1neighmatrix, $lx, $ly, $isHash, $indx, $indy, $ncmatrix, $xmax, $ymax) = @_;

	my $isPeriodicX = $isHash->{"isPeriodicX"};
	my $isPeriodicY = $isHash->{"isPeriodicY"};
	my $isArmchairX = $isHash->{"isArmchairX"};
	my $isLeft = $isHash->{"isLeft"};

	#for NN bonds X=(1,1), Y=(1,-1), Z=(0,1):
	my @dirx = (1, 1, 0);
	my @diry = (1, -1, 1);
	my @ylimit = (1, 0, 1);
	my $n = scalar(@$indx);

	for (my $dir = 0; $dir < 3; ++$dir) {
		for (my $i = 0; $i < $n; ++$i) {
			my $ix = int($indx->[$i]);
			my $iy = int($indy->[$i]);

			my $mx = int($ix + $dirx[$dir]); # 1 1 0
			my $my = int($iy + $diry[$dir]); # 1 -1 1

			if ($mx < $xmax + $dirx[$dir] and $my < $ymax + $ylimit[$dir] and $my >= 0 and defined(getMatrixElem($ncmatrix, $mx, $my))) {
				my $j = int(getMatrixElem($ncmatrix, $mx, $my));
				setMatrixElem($n1neighmatrix, $i, $dir, $j);
				setMatrixElem($n1neighmatrix, $j, $dir, $i);
			}

			if ($dir == 2) {
				if ($isPeriodicY and !$isLeft) {
					$mx = int($ix - $ly);
					$my = 0;
					if ($mx < $xmax and $iy == $ymax and defined(getMatrixElem($ncmatrix, $mx, $my))) {
						my $j = int(getMatrixElem($ncmatrix, $mx, $my));
						setMatrixElem($n1neighmatrix, $i, 2, $j);
						setMatrixElem($n1neighmatrix, $j, 2, $i);
					}
				} elsif ($isPeriodicY and $isLeft) {
					$mx = int($ix + $ly);
					$my = 0;
					if($mx < $xmax + 1 and $iy == $ymax and defined(getMatrixElem($ncmatrix, $mx, $my))) {
						my $j = int(getMatrixElem($ncmatrix, $mx, $my));
						setMatrixElem($n1neighmatrix, $i, 2, $j);
						setMatrixElem($n1neighmatrix, $j, 2, $i);
					}
				}
			} elsif ($dir == 1) {
				if ($isPeriodicX and !$isLeft and !($iy & 1) and getMatrixElem($ncmatrix, $mx, $my) < $ly*2) {
					$mx = int($ix + $lx*2 - 1);
					$my = int($iy + 1);

					if( $mx < $xmax + 1 and $iy < $ymax + 1 and defined(getMatrixElem($ncmatrix, $mx, $my))) {
						my $j = int(getMatrixElem($ncmatrix, $mx, $my));
						setMatrixElem($n1neighmatrix, $i, 1, $j);
						setMatrixElem($n1neighmatrix, $j, 1, $i);
					}
				}
			} elsif ($dir == 0) {
				if ($isPeriodicX and $isLeft and ($iy & 1) and getMatrixElem($ncmatrix, $mx, $my) < $ly*2) {
					$mx = int($ix + $lx*2 - 1);
					$my = int($iy - 1);
					if ($mx < $xmax + 1 and $iy < $ymax + 1 and defined(getMatrixElem($ncmatrix, $mx, $my))) {
						my $j = int(getMatrixElem($ncmatrix, $mx, $my));
						setMatrixElem($n1neighmatrix, $i, 0, $j);
						setMatrixElem($n1neighmatrix, $j, 0, $i);
					}
				}
			}
		}
	}
}


sub honeycombZigzagSetThirdNeigh
{
	my ($n3neighmatrix, $lx, $ly, $isHash, $indx, $indy, $ncmatrix, $xmax, $ymax) = @_;

	my $isPeriodicX = $isHash->{"isPeriodicX"};
	my $isPeriodicY = $isHash->{"isPeriodicY"};
	my $isArmchairX = $isHash->{"isArmchairX"};
	my $isLeft = $isHash->{"isLeft"};
	#for NNNN bonds X3=(-2,-1), Y3=(+2,-1), Z3=(0,-3):
	my @dirx = (2,  2,  0);
	my @diry = (1,  -1,  3);
	my @ylimit = (1, 0, 3);
	my $n = scalar(@$indx);
	
	if($isPeriodicX) {
		print("\n\nPeriodicX not yet supported with J3.\n\n\n");
	}

	for (my $dir = 0; $dir < 3; ++$dir) {
		for (my $i = 0; $i < $n; ++$i) {
			my $ix = int($indx->[$i]);
			my $iy = int($indy->[$i]);

			my $mx = int($ix + $dirx[$dir]); # 2  2 0
			my $my = int($iy + $diry[$dir]); # 1 -1 3

			if ($mx < $xmax + $dirx[$dir] and $my < $ymax + $ylimit[$dir] and $my >= 0 and defined(getMatrixElem($ncmatrix, $mx, $my))) {
				my $j = int(getMatrixElem($ncmatrix, $mx, $my));
				setMatrixElem($n3neighmatrix, $i, $dir, $j);
				setMatrixElem($n3neighmatrix, $j, $dir, $i);
			}

			if ($dir == 2) {
				if ($isPeriodicY and !$isLeft) {
					$mx = int($ix - $ly);
					$my = 1;
					if ($mx < $xmax and $iy == (-1+$ymax) and defined(getMatrixElem($ncmatrix, $mx, $my))) {
						my $j = int(getMatrixElem($ncmatrix, $mx, $my));
						setMatrixElem($n3neighmatrix, $i, 2, $j);
						setMatrixElem($n3neighmatrix, $j, 2, $i);
					}
				} elsif ($isPeriodicY and $isLeft) {
					$mx = int($ix + $ly);
					$my = 1;
					if($mx < $xmax + 1 and $iy == (-1+$ymax) and defined(getMatrixElem($ncmatrix, $mx, $my))) {
						my $j = int(getMatrixElem($ncmatrix, $mx, $my));
						setMatrixElem($n3neighmatrix, $i, 2, $j);
						setMatrixElem($n3neighmatrix, $j, 2, $i);
					}
				}
			} elsif ($dir == 1) {
				if ($isPeriodicY and !$isLeft) {
					$mx = int($ix + $ly + 2);
					$my = int($ymax);
					if( $mx < $xmax + 1 and ($iy == 0) and defined(getMatrixElem($ncmatrix, $mx, $my))) {
						my $j = int(getMatrixElem($ncmatrix, $mx, $my));
						setMatrixElem($n3neighmatrix, $i, 1, $j);
						setMatrixElem($n3neighmatrix, $j, 1, $i);
					}
				}
			} elsif ($dir == 0) {
				if ($isPeriodicY and !$isLeft) {
					$mx = int($ix - $ly + 2);
					$my = int(0);
					if( $mx >= 0 and ($iy == $ymax) and defined(getMatrixElem($ncmatrix, $mx, $my))) {
						my $j = int(getMatrixElem($ncmatrix, $mx, $my));
						setMatrixElem($n3neighmatrix, $i, 0, $j);
						setMatrixElem($n3neighmatrix, $j, 0, $i);
					}
				}
			} 
		}
	}
}


sub honeycombArmchairSetNeigh
{
	my ($n1neighmatrix, $lx, $ly, $isHash, $indx, $indy, $ncmatrix, $xmax, $ymax) = @_;

	my $isPeriodicX = $isHash->{"isPeriodicX"};
	my $isPeriodicY = $isHash->{"isPeriodicY"};
	my $isArmchairX = $isHash->{"isArmchairX"};
	my $isLeft = $isHash->{"isLeft"};

	#for NN bonds X=..., Y=..., Z=...:
	my @dirx = (1, 1, 1);
	my @diry = (-1, 1, 0);
	my @xlimit = (1, 1, 1);
	my @ylimit = (1, 1, 0);
	my $n = scalar(@$indx);
	print "\nxmax=", $xmax;
	print "\nymax=", $ymax;
	
	for (my $dir = 0; $dir < 3; ++$dir) {
		print "\n=====NEW DIRECTION=====\n";
		for (my $i = 0; $i < $n; ++$i) {

			my $ix = int($indx->[$i]);
			my $iy = int($indy->[$i]);
			my $mx = int($ix + $dirx[$dir]); # 
			my $my = int($iy + $diry[$dir]); # 
			print "\n";
			print $i, " ", $ix, " ", $iy, " ", $mx, " ", $my;
			print "\n";

			if ($mx < $xmax + $dirx[$dir]  and $my <= $ymax + $ylimit[$dir] and $my >= 0 and defined(getMatrixElem($ncmatrix, $mx, $my))) {
				my $j = int(getMatrixElem($ncmatrix, $mx, $my));
				print $i, " ", $j, "\n";
				setMatrixElem($n1neighmatrix, $i, $dir, $j);
				setMatrixElem($n1neighmatrix, $j, $dir, $i);
			}
			
			if ($dir == 0) {
				if ($isPeriodicY) {
					$mx = int($ix-1);
					$my = 0;
					if ($mx < $xmax and $iy == ($ymax) and defined(getMatrixElem($ncmatrix, $mx, $my))) {
						my $j = int(getMatrixElem($ncmatrix, $mx, $my));
						print "\nheyo0, ix=",$ix, " iy=",$iy," mx=",$mx," my=",$my,"\n";
						print $i, " ", $j, "\n";
						setMatrixElem($n1neighmatrix, $i, 0, $j);
						setMatrixElem($n1neighmatrix, $j, 0, $i);
					}
				}
			} elsif ($dir == 1) {
				if ($isPeriodicY) {
					$mx = int($ix+1);
					$my = 0;
					if ($mx <= $xmax and $iy == ($ymax) and defined(getMatrixElem($ncmatrix, $mx, $my))) {
						my $j = int(getMatrixElem($ncmatrix, $mx, $my));
						print "\nheyo1, ix=",$ix, " iy=",$iy," mx=",$mx," my=",$my,"\n";
						print $i, " ", $j, "\n";
						setMatrixElem($n1neighmatrix, $i, 1, $j);
						setMatrixElem($n1neighmatrix, $j, 1, $i);
					}
				}
			} elsif ($dir == 2) {
				if ($isPeriodicX and $ix>0) {
					$mx = 0;
					$my = int($iy);
					if ($ix == $xmax and $iy < ($ymax) and defined(getMatrixElem($ncmatrix, $mx, $my))) {
						my $j = int(getMatrixElem($ncmatrix, $mx, $my));
						print "\nheyo1, ix=",$ix, " iy=",$iy," mx=",$mx," my=",$my,"\n";
						print $i, " ", $j, "\n";
						setMatrixElem($n1neighmatrix, $i, 2, $j);
						setMatrixElem($n1neighmatrix, $j, 2, $i);
					}
				}
			}

			
		}
	}
}


sub honeycombArmchairSetThirdNeigh
{
	my ($n3neighmatrix, $lx, $ly, $isHash, $indx, $indy, $ncmatrix, $xmax, $ymax) = @_;

	my $isPeriodicX = $isHash->{"isPeriodicX"};
	my $isPeriodicY = $isHash->{"isPeriodicY"};
	my $isArmchairX = $isHash->{"isArmchairX"};
	my $isLeft = $isHash->{"isLeft"};
	#for NNNN bonds X3=(+1,-2), Y3=(+1,+2), Z3=(+3,0):
	my @dirx = (1,   1,  3);
	my @diry = (-2,  2,  0);
	my @ylimit = (2, 1, 1);
	my $n = scalar(@$indx);	

	for (my $dir = 0; $dir < 3; ++$dir) {
		for (my $i = 0; $i < $n; ++$i) {
			my $ix = int($indx->[$i]);
			my $iy = int($indy->[$i]);

			my $mx = int($ix + $dirx[$dir]); # 1   1  3
			my $my = int($iy + $diry[$dir]); # -2  2  0

			if ($mx < $xmax + $dirx[$dir] and $my < $ymax + $ylimit[$dir] and $my >= 0 and defined(getMatrixElem($ncmatrix, $mx, $my))) {
				my $j = int(getMatrixElem($ncmatrix, $mx, $my));
				setMatrixElem($n3neighmatrix, $i, $dir, $j);
				setMatrixElem($n3neighmatrix, $j, $dir, $i);
			}

			if ($dir == 2) {
				if ($isPeriodicX and ($ix == 1 ) ){
					$mx = int(4*$lx - 2);
					$my = $iy;
					my $j = int(getMatrixElem($ncmatrix, $mx, $my));
					setMatrixElem($n3neighmatrix, $i, 2, $j);
					setMatrixElem($n3neighmatrix, $j, 2, $i);
				}
			} elsif ($dir == 1) {
				if ($isPeriodicX and $ix == $xmax and $my < $ymax + $ylimit[$dir] and $my >= 0) {
					$mx = 0;
					my $j = int(getMatrixElem($ncmatrix, $mx, $my));
					setMatrixElem($n3neighmatrix, $i, 1, $j);
					setMatrixElem($n3neighmatrix, $j, 1, $i);
				} elsif ($isPeriodicY and $iy == $ymax and ($ix & 1) and $mx < $xmax and $mx >= 0){
					$my = 1;
					my $j = int(getMatrixElem($ncmatrix, $mx, $my));
					setMatrixElem($n3neighmatrix, $i, 1, $j);
					setMatrixElem($n3neighmatrix, $j, 1, $i);
				} elsif ($isPeriodicY and $iy == $ymax-1 and ($ix & 1) and $mx < $xmax and $mx > 0){
					$my = 0;
					my $j = int(getMatrixElem($ncmatrix, $mx, $my));
					setMatrixElem($n3neighmatrix, $i, 1, $j);
					setMatrixElem($n3neighmatrix, $j, 1, $i);
				} elsif ($isPeriodicX and $isPeriodicY and $iy >= $ymax-1 and ($ix & 1) and $ix >= $xmax -1){
					$mx = int($ix + 1 - 4*$lx);
					$my = int($iy + 2 - 2*$ly);
					my $j = int(getMatrixElem($ncmatrix, $mx, $my));
					setMatrixElem($n3neighmatrix, $i, 1, $j);
					setMatrixElem($n3neighmatrix, $j, 1, $i);
				}
				
			} elsif ($dir == 0) {
				if ($isPeriodicX and $ix == $xmax and $my < $ymax + $ylimit[$dir] and $my >= 0) {
					$mx = 0;
					my $j = int(getMatrixElem($ncmatrix, $mx, $my));
					setMatrixElem($n3neighmatrix, $i, 0, $j);
					setMatrixElem($n3neighmatrix, $j, 0, $i);
				} elsif ($isPeriodicY and $iy >= $ymax-1 and !($ix & 1) and $ix < $xmax and $ix-1 >= 0){
					$mx = int($ix -1);
					$my = int($iy + 2 - 2*$ly);
					my $j = int(getMatrixElem($ncmatrix, $mx, $my));
					setMatrixElem($n3neighmatrix, $i, 0, $j);
					setMatrixElem($n3neighmatrix, $j, 0, $i);
				} elsif ($isPeriodicX and $isPeriodicY and $iy-2 < 0){
					$mx = int($ix + 1 - 4*$lx);
					$my = int($iy - 2 + 2*$ly);
					my $j = int(getMatrixElem($ncmatrix, $mx, $my));
					setMatrixElem($n3neighmatrix, $i, 0, $j);
					setMatrixElem($n3neighmatrix, $j, 0, $i);
				}
			}
		}
	}
}

sub honeycombZigzagSetSiteCoordinates
{
	my ($tindx, $tindy, $ptsx, $ptsy, $ncmatrix, $lx, $ly, $isHash) = @_;

	my $isPeriodicX = $isHash->{"isPeriodicX"};
	my $isPeriodicY = $isHash->{"isPeriodicY"};
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

	my $isPeriodicX = $isHash->{"isPeriodicX"};
	my $isPeriodicY = $isHash->{"isPeriodicY"};
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


sub honeycombZigzagPlot
{
	my ($indx, $indy, $tindx, $tindy, $n1neighmatrix, $ncmatrix, $isHash) = @_;
	my $str = "";
	my $tindyMax = max @$tindy;
	my $tindyMin = min @$tindy;
	
	for (my $i = 0; $i < $n; ++$i) {
		
		my $sqix = int($indx->[$i]);
		my $sqiy = int($indy->[$i]);
		my $ix = $tindx->[$i];
		my $iy = $tindy->[$i];
		my $ymove = 1.0;
	
		my $e0 = getMatrixElem($n1neighmatrix, $i, 0);
		if (defined($e0)) {
			my $jx = $tindx->[$e0];
			my $jy = $tindy->[$e0];
			my $enc = getMatrixElem($ncmatrix, $sqix, $sqiy);
			if ($isHash->{"isLeft"}) {
				if ($enc < $ly*2 or $enc > $n - $ly*2 - 1) {
					$str .= plotLine($ix, $jx, $iy, $jy, "mystyle1");
				} else { 
					$str .= plotLine($ix, $jx, $iy, $jy, "mystyle2");
				}
			} else {
				$str .= plotLine($ix, $jx, $iy, $jy, "mystyle2");
			}
		}

		my $e1 = getMatrixElem($n1neighmatrix, $i, 1);
		if( defined($e1)) {
			my $jx = $tindx->[$e1];
			my $jy = $tindy->[$e1];

			if (!$isHash->{"isLeft"}) {
				if($i < $ly*2 - 1 and !($i & 1)) {
					$str .= plotLine($ix, $jx, $iy, $jy, "mystyle3");
				} elsif ($i < $ly*2 - 1 or $i > $n - $ly*2) {
					$str .= plotLine($ix, $jx, $iy, $jy, "mystyle3");
				} else { 
					$str .= plotLine($ix, $jx, $iy, $jy, "mystyle4");
				}
			} else {
				$str .= plotLine($ix, $jx, $iy, $jy, "mystyle4");
			}
		}

		my $e2 = getMatrixElem($n1neighmatrix, $i, 2);
		if( defined($e2)) {
			my $jx = $tindx->[$e2];
			my $jy = $tindy->[$e2];
			if ($iy == $tindyMax) {
				$str .= plotLine($ix, $jx, $iy, $jy, "mystyle5") if ($dashed);
			} elsif ($iy == $tindyMin) {
				$str .= plotLine($ix, $jx, $iy, $jy, "mystyle5") if ($dashed);
			} else { 
				$str .= plotLine($ix, $jx, $iy, $jy, "mystyle6");
			}
		}
	}

	for (my $i = 0; $i < $n; ++$i) {
		
		my $sqix = int($indx->[$i]);
		my $sqiy = int($indy->[$i]);
		my $ix = $tindx->[$i];
		my $iy = $tindy->[$i];

		$str .= plotNode($ix, $iy, 0.07, "$i");
	}

	return $str;
}


sub honeycombArmchairPlot
{
	my ($indx, $indy, $tindx, $tindy, $n1neighmatrix, $ncmatrix, $isHash) = @_;
	my $str = "";
	my $tindyMax = max @$tindy;
	my $tindyMin = min @$tindy;
	my $tindxMax = max @$tindx;
	my $tindxMin = min @$tindx;
	my $isArmchairX = $isHash->{"isArmchairX"};
	if ($isArmchairX){
		$n = 4*$lx*$ly;
	}
	
	for (my $i = 0; $i < $n; ++$i) {
		
		my $sqix = int($indx->[$i]);
		my $sqiy = int($indy->[$i]);
		my $ix = $tindx->[$i];
		my $iy = $tindy->[$i];
		my $ymove = 1.0;
	
		my $e0 = getMatrixElem($n1neighmatrix, $i, 0);
		if (defined($e0)) {
			my $jx = $tindx->[$e0];
			my $jy = $tindy->[$e0];
			my $enc = getMatrixElem($ncmatrix, $sqix, $sqiy);
			if (($iy == $tindyMax and ($i & 1)) or ($iy == $tindyMin and !($i & 1))){
				$str .= plotLine($ix, $jx, $iy, $jy, "mystyle1");
			} else {
				$str .= plotLine($ix, $jx, $iy, $jy, "mystyle2");
			}
		}

		my $e1 = getMatrixElem($n1neighmatrix, $i, 1);
		if( defined($e1)) {
			my $jx = $tindx->[$e1];
			my $jy = $tindy->[$e1];
			if (($iy == $tindyMax and !($i & 1)) or ($iy == $tindyMin and ($i & 1))){
				$str .= plotLine($ix, $jx, $iy, $jy, "mystyle3");
			} else {
				$str .= plotLine($ix, $jx, $iy, $jy, "mystyle4");
			}
		}

		my $e2 = getMatrixElem($n1neighmatrix, $i, 2);
		if( defined($e2)) {
			my $jx = $tindx->[$e2];
			my $jy = $tindy->[$e2];
			print "\ni=",$i,", jx=",$jx,", jy=",$jy;
			if ($ix == $tindxMax or $ix == $tindxMin) {
				$str .= plotLine($ix, $jx, $iy, $jy, "mystyle7");
			} else { 
				$str .= plotLine($ix, $jx, $iy, $jy, "mystyle6");
			}
		}
	}

	for (my $i = 0; $i < $n; ++$i) {
		
		my $sqix = int($indx->[$i]);
		my $sqiy = int($indy->[$i]);
		my $ix = $tindx->[$i];
		my $iy = $tindy->[$i];

		$str .= plotNode($ix, $iy, 0.07, "$i");
	}

	return $str;
}

sub plotLine
{
	my ($ix, $jx, $iy, $jy, $style) = @_;
	my $cix = sprintf("%.2f", $ix);
	my $ciy = sprintf("%.2f", $iy);
	my $cjx = sprintf("%.2f", $jx);
	my $cjy = sprintf("%.2f", $jy);
	return "\\draw[$style] ($cix, $ciy) -- ($cjx, $cjy);\n"
}

sub plotNode
{
	my ($ix, $iy, $radius, $label) = @_;
	my $cix = sprintf("%.2f", $ix);
	my $ciy = sprintf("%.2f", $iy);
	my $prespx = ($label < 10) ? 0.2 : 0.3;
	my $spx = ($label & 1) ?  $prespx : -0.28;
	my $spy = ($label & 1) ? 0.02 : -0.02;
	my $str = "\\draw[mystycir] ($cix, $ciy) circle ($radius);\n";
	if ($withNumbers) {
		$str .= "\\node at ($cix + $spx, $ciy + $spy) {\\tiny $label };\n";
	}

	return $str;
}

sub getMatrixElem
{
	my ($a, $mx, $my) = @_;
	return $a->{"data"}->[$mx + $my*$a->{"rows"}];
}

sub setMatrixElem
{
	my ($a, $mx, $my, $val) = @_;
	$a->{"data"}->[$mx + $my*$a->{"rows"}] = $val;
}

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

