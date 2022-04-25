#!/usr/bin/perl

use strict;
use warnings;
use Math::Trig;
use lib ".";
use OmegaUtils;
use Getopt::Long qw(:config no_ignore_case);

my $usage = "-f dollarizedInput [-M mForQ] [-S site] [-p] [-r] [-z]\n";

my ($templateInput,$site,$m,$GlobalNumberOfSites);
my ($siteForSpectrum,$mForQ,$mMax,$wantsRealPart);
my $zeroAtCenter = 0;
my $isPeriodic = 0;
my $nonNegativeOnly = 0;
my $noFourier = 0;
my $comment = "#";

GetOptions('f=s' => \$templateInput,
           'S:i' => \$siteForSpectrum,
           'm:i' => \$mForQ,
           'p' => \$isPeriodic,
           'M:i' => \$mMax,
		   'c=s' => \$comment,
           'r' => \$wantsRealPart,
		   'N' => \$nonNegativeOnly,
           'z' => \$zeroAtCenter,
           'X' => \$noFourier) or die "$usage\n";

(defined($templateInput) and defined($isPeriodic)) or die "$0: USAGE: $usage\n";

my ($omega0, $omegaTotal, $omegaStep, $centralSite);
my $geometryName;
my $geometrySubName = "NONE";
my $geometryLeg = 1;
my $orbitals = 1;
my $omegaOffset = 0;
my $jacksOrLorentz = "none";
my $ChebyC = 0;
my $ChebySign = 1;
my $testoutputfile = "runForinput\$";
my ($lx, $ly) = (0, 0);
my $options = "";

my $isAinur = OmegaUtils::isAinur  ($templateInput);

my $hptr = {"${comment}OmegaBegin" => \$omega0,
            "${comment}OmegaTotal" => \$omegaTotal,
            "${comment}OmegaStep" => \$omegaStep,
            "${comment}OmegaOffset" => \$omegaOffset,
            "#Lx" => \$lx,
            "#Ly" => \$ly,
            "#options" => \$options,
            "GeometryKind" => \$geometryName,
            "GeometrySubKind" => \$geometrySubName,
            "LadderLeg" => \$geometryLeg,
            "Orbitals" => \$orbitals,
            "${comment}ChebyshevC" => \$ChebyC,
            "${comment}ChebyshevSign" => \$ChebySign,
            "${comment}JacksonOrLorentz" => \$jacksOrLorentz,
            "TotalNumberOfSites" => \$GlobalNumberOfSites,
            "OutputFile" => \$testoutputfile};

OmegaUtils::getLabels($hptr, $templateInput);

$centralSite = getCentralSite($templateInput, $isAinur);

if ($isAinur) {
	$geometryName =~ s/[\";]//g;
	$GlobalNumberOfSites =~ s/;//g;
	$geometrySubName =~ s/[\";]//g;
	$omega0 =~ s/[\";]//g;
	$omegaTotal =~ s/[\";]//g;
	$omegaStep =~ s/[\";]//g;
	$omegaOffset =~ s/[\";]//g;
	$geometryLeg =~ s/;//;
}

$hptr->{"isPeriodic"} = $isPeriodic;
$hptr->{"mMax"} = $mMax;
$hptr->{"centralSite"} = $centralSite;
$hptr->{"isCheby"} = findIfWeAreCheby($jacksOrLorentz, $testoutputfile, $ChebyC);

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
open(SPACEOUT, ">", "out.space") or die "$0: Cannot write to out.space : $!\n";
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
print STDERR "$0: Momentum Spectrum written to $outSpectrum\n";
close(SPACEOUT);
print STDERR "$0: Spatial Spectrum written to out.space\n";
exit(0) if ($noFourier);
my $wantsRealOrImag = (defined($wantsRealPart)) ? "real" : "imag";
my $omegaMax = $omega0 + $omegaStep * $omegaTotal;

OmegaUtils::printGnuplot($outSpectrum, $geometry, $isPeriodic, $zeroAtCenter, $nonNegativeOnly);

close(LOGFILEOUT);
print STDERR "$0: Log written to $logFile\n";

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

	writeSpaceValues($omega, \@spaceValues);

	return if ($noFourier);

	my @qValues;
	OmegaUtils::fourier(\@qValues,\@spaceValues,$geometry,$hptr);
	print LOGFILEOUT "$0: Number of k values ".scalar(@qValues)."\n";
	OmegaUtils::writeFourier($array,\@qValues,$geometry);
}

sub writeSpaceValues
{
	my ($omega, $array) = @_;
	my $n = scalar(@$array);
	print SPACEOUT "$omega $n\n";
	for (my $i = 0; $i < $n; ++$i) {
		my $ptr = $array->[$i];
		my $vv1 = $ptr->[0];
		my $vv2 = $ptr->[1];
		print SPACEOUT "$i $vv1 $vv2\n";
	}
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
		my $thisLine = $_;
		foreach my $label (@$labels) {
			my $isGs = ($thisLine =~ /gs/ or $thisLine =~ /X0/);
			next unless (/$label/ and $isGs);
			$status = $label;
			$skip = 0;
		}

		next if ($skip);

		chomp;
		my @temp = split;
		die "$0: Line $_\n" unless (scalar(@temp)==5);

		my $site = $temp[0];
		my $time = ($hptr->{"isCheby"}) ? $temp[2] : 0;
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
		if (!defined($vv2)) {
			print STDERR "$0: P2: Undefined value for site = $i and omega = $omega\n";
			$vv2 = 0.0;
		}

		if (isComplexNumber($vv2)) {
			($vv1, $vv2) = complexNumberToList($vv2);
		}
		
		if (!defined($vv1)) {
			print STDERR "$0: P3: Undefined value for site = $i and omega = $omega\n";
			$vv1 = 0.0;
		}
		
		$array->[$i] = [$vv1, $vv2];
	}
}

sub isComplexNumber
{
	$_ = shift;
	if (/^ *\(/) {
		if (/\) *$/) {
			return 1;
		}
		
		die "$0: $_ isn't real or complex!?\n";
	}
	
	return 0;
}

sub complexNumberToList
{
	my ($x) = @_;
	$_ = $x;
	if (s/^ *\(//) {
		if (s/\) *$//) {
			my @temp = split/,/;
			scalar(@temp) == 2 or notAComplexNumber($x);
			return @temp;
		}
	}
	
	die "$0: $x isn't real or complex!?\n";
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
	my ($dir, $base) = getDirAndBase($file);
	$file = "runFor$base.cout";

	open(FIN, "<", "$file") or die "$0: Cannot open $file : $!\n";

	correctionVectorReadOpen([$v1], ["P1"], $file, \*FIN);
	close(FIN);
}

sub getDirAndBase
{
	my ($file) = @_;
	my $b1 = ($file =~ /^\//);
	my $b2 = ($file =~ /^\.\.\//);
	return ("./", $file) if (!$b1 && !$b2);

	my @temp = split/\//, $file;
	my $n = scalar(@temp);
	die "$0: getDirAndBase: expected n >= 2, got $n\n" if ($n < 2);

	my $dir = $temp[0]."/";
	my $base = $temp[$n - 1];
	for (my $i = 1; $i < $n - 1; ++$i) {
		$dir = $temp[$i]."/";
	}

	return ($dir, $base);
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
	my $epsilont = 0.025;
	my $Wprime = 1.0-0.5*$epsilont;
	my (@factor, @polyXfactor);

	for (my $om = 0; $om < $omegaTotal; ++$om) {
		my $omega = $omega0+$om*$omegaStep;
		# Scaling
		my $omegaPrime = ($ChebySign*$omega-$omega0)*$ChebyC-$Wprime;
		my $underSqrt = 1.0-$omegaPrime*$omegaPrime;
		die "$0: $omegaPrime is such that 1-w'^2 < 0\n" if ($underSqrt < 0);
		$factor[$om] = 1.0/(pi*sqrt($underSqrt));
		for (my $n = 1; $n < $times; ++$n) {
			$polyXfactor[$n + $times*$om] = $ChebyC*$factor[$om]*cos($n*acos($omegaPrime));
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

	return if ($noFourier);
	my @qValues;
	OmegaUtils::fourier(\@qValues,\@spaceValues,$geometry,$hptr);

	my  @array;
	OmegaUtils::writeFourier(\@array,\@qValues,$geometry);
	die "doCheby: array is empty\n" if (scalar(@array) == 0);
	printSpectrum(\@array);
}

sub procAllQs
{
	my ($ind, $omega, $centralSite, $geometry) = @_;
	my @array;
	procCommon(\@array, $ind, $omega, $centralSite, $geometry);
	return if ($noFourier);
	die "procAllQs: array is empty\n" if (scalar(@array) == 0);
	printSpectrum(\@array);
}

sub execThis
{
	my ($cmd) = @_;
	print LOGFILEOUT "$0: About to execute $cmd\n";
	system($cmd);
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

sub findIfWeAreCheby
{
	my ($jacksOrLorentz, $testoutputfile, $chebyC) = @_;
	my $b1 = ($jacksOrLorentz eq "none");
	my $b2 = ($testoutputfile =~ /\$/);

	return 0 if ($b1 and $b2);

	if ($chebyC == 0) {
		die "$0: Needs ChebyshevC in input\n";
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

sub getCentralSite
{
	my ($input, $isAinur) = @_;
	open(FILE, "<", $input) or die "$0: Cannot open $input : $!\n";
	my $c;
	while (<FILE>) {
		chomp;
		if (!$isAinur && /^[ \t]*TSPSites[ \t]+1[ \t]+(\d+)[ \t]*$/) {
			$c = $1;
			last;
		}

		if ($isAinur && /^[ \t]*TSPSites[ \t]*=[ \t]*\[[ \t]*(\d+)[ \t]*\]/) {
			$c = $1;
			last;
		}
	}

	close(FILE);
	defined($c) or die "$0: Could not find TSPSites in $input\n";
	return $c;
}

