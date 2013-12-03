#!/usr/bin/perl

use strict;
use warnings;
use Utils;

my ($site,$site2,$what,$root,$cOrN) = @ARGV;
defined($cOrN) or die "USAGE: $0 site site2 what root cOrN\n";
Utils::checkRange($cOrN,"c","n");

my %params = Utils::loadParams("params.txt");

my $templateInput = "inputTemplate.inp";

my $n = Utils::getLabel($templateInput,"TotalNumberOfSites=");
my @spectral;

my $input = createInput($site,$site2);

my $siteMin = ($site < $site2) ? $site : $site2;
my $siteMax = ($site < $site2) ? $site2 : $site;
my $output = "$root${siteMin}_$siteMax";

my $b = ($what & 2);
if ($b && ($site2 >= $site)) {
	print STDERR "Running ./lanczos -f $input  -g $cOrN &> $output.comb\n";
	system("./lanczos -f $input  -g $cOrN &> $output.comb");
}

die "$0: $output.comb does not exist\n" unless (-r "$output.comb");

$b = ($what & 1);

if ($b  && ($site2 >= $site)) {
	my ($psimagLite,$begin,$end,$step,$eps) = ($params{"PsimagLite"},$params{"OmegaBegin"},$params{"OmegaEnd"},
		                                   $params{"OmegaStep"},$params{"OmegaEps"});
	system("$psimagLite/drivers/continuedFractionCollection -f $output.comb -b $begin -e $end -s $step -d $eps > $output.cf");
}

die "$0: $output.cf does not exist\n" unless (-r "$output.cf");

print STDERR "$0: Done $site $site2\n";

sub createInput
{
	my @sites=@_;
	my $type = 0;
	my $file="input$type.inp";
	open(FOUT,">$file") or die "$0: Cannot write to $file\n";

	my @matrix;
	for (my $i = 0; $i < 4; $i++) {
		$matrix[$i] = zeroMatrix();
	}

	$matrix[0] = findMatrix($type,0);
	$matrix[1] = findMatrix($type,1);

	my $sites = "2 @sites";
	my $loops = "2 0 0";
	my $U = Utils::getLabel($templateInput,"##U=");
	my $hubbardU = setVector($n,$U);
	my $V = Utils::getLabel($templateInput,"##V=");
	my $potentialV = setVector(2*$n,$V);
	my $data = "data.txt";

	open(FILE,"$templateInput") or die "$0: Cannot open $templateInput: $!\n";

	while(<FILE>) {
		next if (/^#/);
		if (/\$([a-zA-Z0-9\[\]]+)/) {
				my $name = $1;
				my $str = "\$".$name;
				my $val = eval "$str";
				defined($val) or die "$0: Undefined substitution for $name\n";
				s/\$\Q$name/$val/g;
		}
		print FOUT;
	}

	close(FILE);
	close(FOUT);

	return $file;
}

sub findMatrix
{
	my ($type,$type2)=@_;
	my $matrix = "0 1 0 0\n0 0 0 0\n0 0 0 1\n0 0 0 0\n";
	if ($type & 1) {
		 $matrix = "0 0 0 0\n1 0 0 0\n0 0 0 0\n0 0 1 0\n";
	}
	return $matrix if ($type2==0);
	if ($type>1) {
		$matrix = multiplyMatrix($matrix,-1);
	}
	return $matrix;
}

sub multiplyMatrix
{
	my ($matrix)=@_;
	$matrix=~s/\n/ /g;
	my @temp=split/ /,$matrix;
	my $matrix2="";
	for (my $i=0;$i<scalar(@temp);$i++) {
		my $_ = -$temp[$i];
		$matrix2 .= $_." ";
	}
	return $matrix2;
}

sub setVector
{
	my ($sites,$U)=@_;
	my $tmp = "$sites ";
	for (my $i=0;$i<$sites;$i++) {
		$tmp .= " $U ";
	}
	return $tmp;
}

sub zeroMatrix
{
	my $id =  "0 0 0 0\n0 0 0 0\n0 0 0 0\n0 0 0 0\n";
	my $matrixEnd = "\nFERMIONSIGN=1\nJMVALUES 0 0\nAngularFactor=1\n\n";
        return "$id$matrixEnd";
}
