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
printGnuplot($outSpectrum, $geometry);

close(LOGFILEOUT);
print STDERR "$0: Log written to $logFile\n";

sub printGnuplot
{
	my ($inFile, $geometry) = @_;
	my $hasPrinted = 0;
	open(FIN, "<", "$inFile") or die "$0: Cannot open $inFile : $!\n";

	my %h;
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
		$h{$omega} = \@temp;

	}

	close(FIN);

	OmegaUtils::printGnuplot(\%h, $geometry, $isPeriodic, $zeroAtCenter);
}

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
	my ($ind, $omega, $centralSite, $geometry) = @_;
	my $n = $GlobalNumberOfSites;
	my $inputRoot = "input";
	my $prefix = "runFor$inputRoot$ind";
	my $outFile = "$prefix.space";
	my $inFile = "$prefix.cout";
	my @values;
	my @values2;

	my $maxSite = correctionVectorRead(\@values,\@values2,$inFile);

	correctionVectorWrite($outFile,\@values,\@values2,$maxSite,$omega);

	$inFile = "$prefix.space";
	$outFile = "$prefix.sq";

	my @spaceValues;
	readSpace(\@spaceValues,$inFile);
	my @qValues;
	OmegaUtils::fourier(\@qValues,\@spaceValues,$geometry,$hptr);
	printFourier($outFile,\@qValues,$geometry);
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
	my ($outFile, $v1, $v2, $maxSite, $omega) = @_;

	open(FOUT, ">", "$outFile") or die "$0: Cannot write to $outFile : $!\n";

	print FOUT "#omega=$omega\n";
	for (my $i = 0; $i < $maxSite; ++$i) {
		my $vv1 = $v1->[$i];
		my $vv2 = $v2->[$i];
		if (!defined($vv1)) {
			print STDERR "$0: Undefined value for site = $i and omega = $omega\n";
			$vv1 = $vv2 = 0.0;
		}
		print FOUT "$i $vv1 $vv2\n";
	}

	close(FOUT);
}

sub readSpace
{
	my ($space,$inFile) = @_;
	my $counter = 0;

	open(FIN, "<", "$inFile") or die "$0: Cannot open $inFile : $!\n";
	while(<FIN>) {
		if (/^#/) {
		        next;
		}

		my @temp=split;
		next unless (scalar(@temp) == 3);
		my @temp2 = ($temp[1],$temp[2]);
		$space->[$temp[0]] = \@temp2;
		die "$0: Line $_\n" unless ($counter == $temp[0]);
		$counter++;
	}

	close(FIN);
	print LOGFILEOUT "$0: Read $counter sites\n";
}

sub procThisOmegaKspace
{
	my ($ind, $omega, $centralSite, $mForQ, $geometry) = @_;
	procCommon($ind, $omega, $centralSite, $geometry);
	my $inputRoot = "input";
	my $prefix = "runFor$inputRoot$ind";
	my $inFile = "$prefix.sq";
	extractValue($inFile,$mForQ);
}

sub procThisOmegaSpace
{
	my ($ind, $omega, $centralSite, $siteForSpectrum) = @_;
	procCommon($ind, $omega, $centralSite, $geometry);
	my $inputRoot = "input";
	my $prefix = "runFor$inputRoot$ind";
	my $inFile = "$prefix.space";
	extractValue($inFile,$siteForSpectrum);
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

	my $outFile = "runForinput$ind.sq";
	printFourier($outFile,\@qValues,$geometry);

	my @array;
	readAllQs(\@array,$ind);
	die "doCheby: array is empty\n" if (scalar(@array) == 0);
	printSpectrum(\@array);
}

sub procAllQs
{
	my ($ind, $omega, $centralSite, $geometry) = @_;
	procCommon($ind, $omega, $centralSite, $geometry);
	my @array;
	readAllQs(\@array,$ind);
	die "procAllQs: array is empty\n" if (scalar(@array) == 0);
	printSpectrum(\@array);
}

sub readAllQs
{
	my ($array,$ind) = @_;
	my $counter = 0;
	my $inputRoot = "input";
	my $prefix = "runFor$inputRoot$ind";
	open(FILE, "<", "$prefix.sq") or die "$0: Cannot open file : $!\n";
	while (<FILE>) {
		chomp;
		my @temp = split;
		my $n = scalar(@temp);
		($n == 3) or next;
		$array->[$counter++] = \@temp;
	}

	close(FILE);
}

sub execThis
{
	my ($cmd) = @_;
	print LOGFILEOUT "$0: About to execute $cmd\n";
	system($cmd);
}

sub printFourier
{
	my ($outFile, $f, $geometry) = @_;
	my $subname = $geometry->{"subname"};

	if ($geometry->{"name"} eq "chain") {
		return printFourierChain($outFile,$f);
	}

	if ($geometry->{"name"} eq "ladder" || $subname eq "average") {
		return printFourierLadder($outFile, $f);
	}

	die "$0: printFourier: undefined geometry ".$geometry->{"name"}."\n";
}

sub printFourierChain
{
	my ($outFile,$f) = @_;

	open(FOUT, ">", "$outFile") or die "$0: Cannot write to $outFile : $!\n";

	my $n = scalar(@$f);
	for (my $m = 0; $m < $n; ++$m) {
		my $q = OmegaUtils::getQ($m, $n, $isPeriodic);
		my $ptr = $f->[$m];
		my @temp = @$ptr;
		print FOUT "$q $temp[0] $temp[1]\n";
	}

	close(FOUT);
}

sub printFourierLadder
{
	my ($outFile, $f) = @_;

	open(FOUT, ">", "$outFile") or die "$0: Cannot write to $outFile : $!\n";

	my $n = scalar(@$f);
	for (my $m = 0; $m < $n; ++$m) {
		my $ptr = $f->[$m];
		my @temp = @$ptr;
		print FOUT "$m @temp\n";
	}

	close(FOUT);
}

sub extractValue
{
	my ($file,$q) = @_;
	open(FILE, "<", $file) or die "$0: Cannot open file $file : $!\n";

	my $omega;
	while(<FILE>) {

		chomp;
		if (/^#omega=(.*$)/) {
			$omega = $1;
			next;
		}

		my @temp = split;
		next unless (scalar(@temp) > 1);
		next unless (abs($temp[0]-$q)<1e-3);
		die "$0: File $file line $_\n" if (!defined($omega));
		print "$omega $temp[1] ";
		print "$temp[2]" if (scalar(@temp) == 3);
		print "\n";
		last;
	}

	close(FILE);
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

