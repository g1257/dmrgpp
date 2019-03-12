#!/usr/bin/perl

use strict;
use warnings;
use Math::Trig;
use lib ".";
use OmegaUtils;
use Getopt::Long qw(:config no_ignore_case);

my $usage = "-f dollarizedInput [-M mForQ] [-S site] [-p] [-r] [-z]\n";

my ($templateInput,$site,$m,$GlobalNumberOfSites);
my ($siteForSpectrum,$mForQ,$isPeriodic,$mMax,$wantsRealPart);
my $zeroAtCenter = 0;

GetOptions('f=s' => \$templateInput,
           'S:i' => \$siteForSpectrum,
           'm:i' => \$mForQ,
           'p' => \$isPeriodic,
           'M:i' => \$mMax,
           'r' => \$wantsRealPart,
           'z' => \$zeroAtCenter) or die "$usage\n";

(defined($templateInput) and defined($isPeriodic)) or die "$0: USAGE: $usage\n";

my ($omega0, $omegaTotal, $omegaStep, $centralSite);
my $geometryName;
my $geometrySubName = "NONE";
my $geometryLeg = 1;
my $orbitals = 1;
my $omegaOffset = 0;
my $jacksOrLorentz = "none";
my $Wstar = 0;
my $epsilont = 0;
my $testoutputfile;

my $hptr = {"#OmegaBegin" => \$omega0,
            "#OmegaTotal" => \$omegaTotal,
            "#OmegaStep" => \$omegaStep,
            "#OmegaOffset" => \$omegaOffset,
            "GeometryKind" => \$geometryName,
            "GeometrySubKind" => \$geometrySubName,
            "LadderLeg" => \$geometryLeg,
            "Orbitals" => \$orbitals,
            "TSPSites 1" => \$centralSite,
            "#JacksonOrLorentz" => \$jacksOrLorentz,
            "TotalNumberOfSites" => \$GlobalNumberOfSites,
            "ChebyshevWstar" => \$Wstar,
            "ChebyshevEpsilon" => \$epsilont,
            "OutputFile" => \$testoutputfile};

OmegaUtils::getLabels($hptr,$templateInput);
$hptr->{"isPeriodic"} = $isPeriodic;
$hptr->{"mMax"} = $mMax;
$hptr->{"centralSite"} = $centralSite;
$hptr->{"isCheby"} = findIfWeAreCheby($jacksOrLorentz, $testoutputfile, $Wstar, $epsilont);

my $logFile = "Log$templateInput";
$logFile =~ s/\..*$//;
$logFile .= ".log";
open(LOGFILEOUT, ">", "$logFile") or die "$0: Cannot write to $logFile : $!\n";

if ($omegaStep < 0) {
	my $beta = -$omegaStep;
	print LOGFILEOUT "$0: Matsubara freq. assumed with beta= $beta\n";
	$omega0 = $omegaStep = 2.0*pi/$beta;
}

my $geometry = {"name" => $geometryName, "leg" => $geometryLeg, "subname" => $geometrySubName};

my $outSpectrum = "out.spectrum";
open(FOUTSPECTRUM, ">", "$outSpectrum") or die "$0: Cannot write to $outSpectrum : $!\n";

my %cheby;
prepareCheby(\%cheby) if ($hptr->{"isCheby"});


for (my $i = $omegaOffset; $i < $omegaTotal; ++$i) {

	my $omega = $omega0 + $omegaStep * $i;
	print FOUTSPECTRUM "$omega ";
	print LOGFILEOUT "$0: About to proc for omega = $omega\n";

	if ($hptr->{"isCheby"}) {
		doCheby($i, \%cheby, $omega, $centralSite, $geometry);
	} elsif (defined($mForQ)) {
		procThisOmegaKspace($i, $omega, $centralSite, $mForQ, $geometry);
	} elsif (defined($siteForSpectrum)) {
		procThisOmegaSpace($i, $omega, $centralSite, $siteForSpectrum, $geometry);
	} else {
		procAllQs($i, $omega, $centralSite, $geometry);
	}

	print LOGFILEOUT "$0: Finished         omega = $omega\n";
}

close(FOUTSPECTRUM);
print STDERR "$0: Spectrum written to $outSpectrum\n";
my $wantsRealOrImag = (defined($wantsRealPart)) ? "real" : "imag";
my $omegaMax = $omega0 + $omegaStep * $omegaTotal;

printSpectrumToColor($outSpectrum,$wantsRealOrImag,$geometry,$omegaMax);
OmegaUtils::printGnuplot($outSpectrum, $geometry, $isPeriodic, $zeroAtCenter);

close(LOGFILEOUT);
print STDERR "$0: Log written to $logFile\n";

sub printSpectrumToColor
{
	my ($inFile,$what,$geometry,$omegaMax) = @_;
	my ($factor, $fileIndices, $leg) = OmegaUtils::getGeometryDetails($geometry);

	foreach my $fileIndex (@$fileIndices) {
		my $outSpectrum = "outSpectrum$fileIndex.color";
		my @colorData;
		my ($counter,$size) = spectrumToColor(\@colorData,
		                                      $inFile,
		                                      $what,
		                                      $geometry,
		                                      $fileIndex);

		open(FOUTSPECTRUM, ">", "$outSpectrum")
			or die "$0: Cannot write to $outSpectrum : $!\n";
		print FOUTSPECTRUM "$counter $size $omegaMax\n";

		my $rows = scalar(@colorData);
		for (my $i = 0; $i < $rows; ++$i) {
			my @thisRow = @{$colorData[$i]};
			my $cols = scalar(@thisRow);
			for (my $j = 0; $j < $cols; ++$j) {
				my $value = int($thisRow[$j]);
				print FOUTSPECTRUM $value." ";
			}

			print FOUTSPECTRUM "\n";
		}

		close(FOUTSPECTRUM);
		print STDERR "$0: Color spectrum written to $outSpectrum\n";
	}
}

sub printSpectrum
{
	my ($array) = @_;

	for (my $j = 0; $j < scalar(@$array); ++$j) {
		my $array2 = $array->[$j];
		my @array2 = @$array2;
		print FOUTSPECTRUM "$array2[1] $array2[2] ";
	}

	print FOUTSPECTRUM "\n";
}

sub procCommon
{
	my ($array, $ind, $omega, $centralSite, $geometry) = @_;
	my $n = $GlobalNumberOfSites;
	my $inputRoot = "input";
	my $prefix = "runFor$inputRoot$ind";
	my $inFile = "$prefix.cout";
	my @values;
	my @values2;

	my $maxSite = correctionVectorRead(\@values,\@values2,$inFile);

	print STDERR "$0: omega=$omega maxSite=$maxSite\n";
	my @spaceValues;
	correctionVectorWrite(\@spaceValues,\@values,\@values2,$maxSite,$omega);

	my @qValues;
	OmegaUtils::fourier(\@qValues,\@spaceValues,$geometry,$hptr);
	writeFourier($array,\@qValues,$geometry);
}

sub correctionVectorRead
{
	my ($v1,$v2,$inFile,$fh) = @_;
	if (-r "$inFile") {
		open(FIN, "<", "$inFile") or die "$0: Cannot open $inFile : $!\n";
	} else {
		my $input = $inFile;
		$input =~s/runFor//;
		$input =~s/\.cout/\.inp/;
		open(FIN,"./toolboxdmrg -f \"$input\" -a grep -E \"<\" 2>/dev/null |")
		or die "$0: Cannot open pipe : $!\n";
	}

	my $maxSite = 0;

	my @v12 = ($v1, $v2);
	my $labels = ["P3", "P2"]; # ORDER IMPORTANT HERE!
	$maxSite = correctionVectorReadOpen(\@v12, $labels, $inFile,\*FIN);
	close(FIN);

	$maxSite++;

	print LOGFILEOUT "$0: correctionVectorRead maxsite= $maxSite\n";
	return $maxSite;
}

sub correctionVectorReadOpen
{
	my ($vs, $labels, $inFile, $fh) = @_;
	my $status = "clear";
	my $maxSite = 0;

	while (<$fh>) {

		next if (/PsiApp\: +CmdLine/);

		my $skip = 1;
		foreach my $label (@$labels) {
			next unless (/$label/ and /gs/);
			$status = $label;
			$skip = 0;
		}

		next if ($skip);

		chomp;
		my @temp = split;
		die "$0: Line $_\n" unless (scalar(@temp)==5);

		my $site = $temp[0];
		my $time = $temp[2];
		my $c = 0;
		foreach my $label (@$labels) {
			++$c;
			next unless ($status eq $label);
			$vs->[$c - 1]->[$site + $time*$GlobalNumberOfSites] = $temp[1];
		}

		$maxSite = $site if ($maxSite < $site);
		$status = "clear";
	}

	return $maxSite;
}

sub correctionVectorWrite
{
	my ($array, $v1, $v2, $maxSite, $omega) = @_;

	for (my $i = 0; $i < $maxSite; ++$i) {
		my $vv1 = $v1->[$i];
		my $vv2 = $v2->[$i];
		if (!defined($vv1)) {
			print STDERR "$0: Undefined value for site = $i and omega = $omega\n";
			$vv1 = $vv2 = 0.0;
		}

		$array->[$i] = [$vv1, $vv2];
	}
}

sub procThisOmegaKspace
{
	my ($ind, $omega, $centralSite, $mForQ, $geometry) = @_;
	my @array;
	procCommon(\@array, $ind, $omega, $centralSite, $geometry);
	my @temp = extractValue(\@array, $mForQ);
	shift @temp;
	print "$omega @temp\n";
}

sub procThisOmegaSpace
{
	my ($ind, $omega, $centralSite, $siteForSpectrum) = @_;
	my @array;
	procCommon(\@array, $ind, $omega, $centralSite, $geometry);
	my @temp = extractValue(\@array ,$siteForSpectrum);
	shift @temp;
	print "$omega @temp\n";
}

sub readCheby
{
	my ($v1, $input) = @_;
	my $file = $input;
	$file =~s/\.inp//;
	$file = "runFor$file.cout";

	open(FIN, "<", "$file") or die "$0: Cannot open $file : $!\n";

	correctionVectorReadOpen([$v1], ["P1"], $file, \*FIN);
	close(FIN);
}

sub dampCheby
{
	my ($times) = @_;

	my @dampG = (1.0);

	if ($jacksOrLorentz eq "Jackson") {
		my $timesPlusOne = ($times + 1.0);
		my $cot = cot(pi/$timesPlusOne);
		for (my $n = 1; $n < $times; ++$n) {
			my $num = ($timesPlusOne - $n) * cos(pi*$n/$timesPlusOne) + sin(pi*$n/$timesPlusOne) * $cot;
			$dampG[$n] = $num/$timesPlusOne;
		}
	} else {
		my $den = sinh(4.0);
		for (my $n = 1; $n < $times; ++$n) {	
			my $num = sinh(4.0*(1.0-$n/$times));
			$dampG[$n] = $num/$den;
		}
	}

	#print STDERR "Damping Factor\n";
	#for (my $i =0; $i < $times; ++$i) {
		#print STDERR "$i\t$dampG[$i]\n";
	#}

	return @dampG;
}

# Does only one omegas
sub chebyRealSpace
{
	my ($om, $cheby) = @_;
	my $factor = $cheby->{"factor"};
	my $ChebyPolyXfactor = $cheby->{"polyXfactor"};
	my $dampG = $cheby->{"dampG"};
	my $vM = $cheby->{"data"};
	my @vMdamped;
	my $times = scalar(@$dampG);

	$vMdamped[0] = [0, 0];
	for (my $p = 1; $p < $GlobalNumberOfSites - 1; ++$p) {
		# adding the zeroth Cheby momentum $n==0
		my $sum = $factor->[$om]*$vM->[$p];

		# adding all the others		
		for (my $n = 1; $n < $times; ++$n) {
		 	$sum += 2.0*$dampG->[$n]*$ChebyPolyXfactor->[$n+$times*$om]*
				                   $vM->[$p + $GlobalNumberOfSites*$n];
		}

		$vMdamped[$p] = [0, $sum]; # 0 would be the real part
	}

	return @vMdamped;
}

sub doFactorAndPolyXfactor
{
	my ($times) = @_;
	my $wPrime = 1.0-0.5*$epsilont;	
	my $scalea = $Wstar/(2.0*$wPrime);
	my (@factor, @polyXfactor);

	for (my $om = 0; $om < $omegaTotal; ++$om) {
        	my $omega = $omega0+$om*$omegaStep;
		# Scaling
		my $omegaPrime = ($omega/$scalea)-$wPrime;
		my $underSqrt = 1.0 - $omegaPrime*$omegaPrime;
		die "$0: $omegaPrime is such that 1-w'^2 < 0 $wPrime $scalea\n" if ($underSqrt < 0); 
		$factor[$om] = 1.0/(pi*sqrt($underSqrt));
		for (my $n = 1; $n < $times; ++$n) {
			$polyXfactor[$n + $times*$om] = $factor[$om]*cos($n*acos($omegaPrime));
		}
	}

	return (\@factor, \@polyXfactor);
}

sub prepareCheby
{
	my ($cheby) = @_;
	my @v; # index on site and time
	readCheby(\@v, $templateInput);
	$cheby->{"data"} = \@v;

	my $times = scalar(@v)/$GlobalNumberOfSites;
	print STDERR "$0: Read times=$times\n";
	my @dampG = dampCheby($times);
	$cheby->{"dampG"} = \@dampG;

	my ($factor, $polyXfactor) = doFactorAndPolyXfactor($times); # returns 2 references
	$cheby->{"factor"} = $factor;
	$cheby->{"polyXfactor"} = $polyXfactor;
}

sub doCheby
{
	my ($ind, $cheby, $omega, $centralSite, $geometry) = @_;
	my @spaceValues = chebyRealSpace($ind, $cheby);

	my @qValues;
	OmegaUtils::fourier(\@qValues,\@spaceValues,$geometry,$hptr);

	my  @array;
	writeFourier(\@array,\@qValues,$geometry);
	die "doCheby: array is empty\n" if (scalar(@array) == 0);
	printSpectrum(\@array);
}

sub procAllQs
{
	my ($ind, $omega, $centralSite, $geometry) = @_;
	my @array;
	procCommon(\@array, $ind, $omega, $centralSite, $geometry);
	die "procAllQs: array is empty\n" if (scalar(@array) == 0);
	printSpectrum(\@array);
}

sub execThis
{
	my ($cmd) = @_;
	print LOGFILEOUT "$0: About to execute $cmd\n";
	system($cmd);
}

sub writeFourier
{
	my ($array, $f, $geometry) = @_;
	my $subname = $geometry->{"subname"};

	if ($geometry->{"name"} eq "chain") {
		return writeFourierChain($array,$f);
	}

	if ($geometry->{"name"} eq "ladder" || $subname eq "average") {
		return writeFourierLadder($array, $f);
	}

	die "$0: writeFourier: undefined geometry ".$geometry->{"name"}."\n";
}

sub writeFourierChain
{
	my ($array, $f) = @_;

	my $n = scalar(@$f);
	for (my $m = 0; $m < $n; ++$m) {
		my $q = OmegaUtils::getQ($m, $n, $isPeriodic);
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

sub extractValue
{
	my ($array, $q) = @_;
	my $narray = scalar(@$array);

	for (my $i = 0; $i < $narray; ++$i) {
		my $ptr = $array->[$i];
		my @temp = @$ptr;
		next unless (scalar(@temp) > 1);
		next unless (abs($temp[0]-$q)<1e-3);
		return @temp;
	}

	die "$0: No q=$q found in array\n";
}

sub spectrumToColor
{
	my ($data,$file,$realOrImag,$geometry,$qyIndex) = @_;
	my $counter = 0;
	my $size;
	my $finalSize;

	open(FIN, "<", $file) or die "$0: Cannot open file $file : $!\n";
	while (<FIN>) {
		next if (/^#/);
		chomp;
		my @temp = split;
		my $n = scalar(@temp);
		next if ($n < 2);
		if (!defined($size)) {
			$size = $n;
			print STDERR "$0: File $file at line $. has set size to $size\n";
		} else {
			($size == $n) or die "$0: Wrong line $_ (expected size $size)\n";
		}

		my @temp2 = getRealOrImagData(\@temp,$realOrImag,$geometry,$qyIndex);
		$finalSize = scalar(@temp2) if (!defined($finalSize));
		$data->[$counter++] = \@temp2;
	}

	close(FIN);

	die "No data found in $file\n" if ($counter == 0);
	print LOGFILEOUT "$0: Read $counter lines size=$size for $realOrImag from $file\n";

	my ($min,$max) = minMaxData($data);
	print LOGFILEOUT "$0: Data min = $min, max = $max\n";

	scaleData($data,$min,$max);
	return ($counter,$finalSize);
}

sub minMaxData
{
	my ($a) = @_;
	my $rows = scalar(@$a);
	die "minMaxData, rows=0\n" if ($rows == 0);
	my ($min,$max) = ($a->[0]->[1],$a->[0]->[1]);
	for (my $i = 0; $i < $rows; ++$i) {
		my @thisRow = @{$a->[$i]};
		my $cols = scalar(@thisRow);
		for (my $j = 0; $j < $cols; ++$j) {
			my $thisValue = $thisRow[$j];
			#next if ($thisValue<0);
			$min = $thisValue if ($min > $thisValue);
			$max = $thisValue if ($max < $thisValue);
		}
	}

	return ($min,$max);
}

sub scaleData
{
	my ($a,$min,$max) = @_;
	my $afactor = 255/($max-$min);
	my $bfactor = -$afactor*$min;
	my $rows = scalar(@$a);
	for (my $i = 0; $i < $rows; ++$i) {
		my @thisRow = @{$a->[$i]};
		my $cols = scalar(@thisRow);
		for (my $j = 0; $j < $cols; ++$j) {
			my $value = $a->[$i]->[$j];
			#if ($value < 0) {
			#	$a->[$i]->[$j] = 0;
			#	next;
			#}

			$a->[$i]->[$j] = $afactor*$value + $bfactor;
		}
	}
}

sub getRealOrImagData
{
	my ($d,$realOrImag,$geometry,$qyIndex) = @_;
	my @temp;
	my $n = scalar(@$d);
	my $start = 1;
	if ($geometry->{"name"} eq "ladder") {
		my $leg = $geometry->{"leg"};
		$n = int($n/$leg);
		$start += $qyIndex*$n;
		$n *= (1+$qyIndex);
	}

	my $j = 0;
	for (my $i = $start; $i < $n; ++$i) {
		next if ($realOrImag eq "imag" && ($i & 1));
		next if ($realOrImag eq "real" && !($i & 1));
		$temp[$j++] = $d->[$i];
	}

	return @temp;
}

sub findIfWeAreCheby
{
	my ($jacksOrLorentz, $testoutputfile, $Wstar, $epsilont) = @_;
	my $b1 = ($jacksOrLorentz eq "none");
	my $b2 = ($testoutputfile =~ /\$/);

	return 0 if ($b1 and $b2);

	if ($Wstar == 0) {
		die "$0: Needs ChebyshevWstar in input\n";
	}

	if ($epsilont == 0) {
		die "$0: Needs ChebyshevEpsilon in input\n";
	}

	return 1 if (!$b1 and !$b2);

	if ($jacksOrLorentz ne "Jackson" and $jacksOrLorentz ne "Lorentz") {
		die "$0: #JacksOrLorentz= Jackson or Lorentz not $jacksOrLorentz\n";
	}

	if ($b1 and !$b2) {
		die "$0: JacksOrLorentz is none, but OutputFile= is NOT dollarized\n";
	}

	die "$0: JacksOrLorentz is NOT none, but OutputFile= is dollarized\n";
}

