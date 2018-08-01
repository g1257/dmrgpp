#!/usr/bin/perl

use strict;
use warnings;
use utf8;

my ($file, $center, $occurrence) = @ARGV;
($file) or die "USAGE: $0 datafile.hd5 [center] [occurrence]\n";

my $label = "/Def/DensityMatrixEigenvalues/Size";
my $n = readDataSet($file, $label);

print "#There are $n center of orthogonalities\n";

my $nodata = "datafile.hd5 doesn't contain the data,";
my $perhaps = "perhaps you didn't add saveDensityMatrixEigenvalues in SolverOptions";

($n > 0) or die "$0: $nodata $perhaps\n";

($center) or exit(0);

if ($center =~ /(^\d+)\+(\d+$)/) {
	$center = $1 - 1;
}

($center =~ /(^\d+$)/) or die "$0: center must be a number or number+number\n";

die "$0: Center too big, must be smaller than $n\n" if ($center >= $n);

$label = "/Def/DensityMatrixEigenvalues/$center/Size";
$n = readDataSet($file, $label);

print "#Center of orthogonality $center has $n occurrences\n";
($n > 0) or die "$0: $nodata for this center of orthogonality\n";

$occurrence = 0 if (!defined($occurrence) and $n == 1);

($occurrence) or exit(0);

($occurrence =~ /(^\d+$)/) or die "$0: occurrence must be a number\n";

die "$0: occurrence too big, must be smaller than $n\n" if ($occurrence >= $n);

$label = "/Def/DensityMatrixEigenvalues/$center/$occurrence";

my @eigs = readDataSet($file, $label);

printVector(\@eigs);

sub readDataSet
{
	my ($file, $label) = @_;

	if ($label =~ /(^[a-zA-Z0-9\/]+)/) {
		$label = $1;
	} else {
		die "$0: Invalid label $label\n";
	}

	if ($file =~ /(^[a-zA-Z0-9\.\/]+)/) {
		$file = $1;
	} else {
		die "$0: Invalid file $file\n";
	}

	die "$0: File $file not readable\n" unless (-r "$file");

	$ENV{"PATH"} = "";
	$ENV{"ENV"} = "";
	$ENV{"BASH_ENV"} = "";
	my $value = "";
	open (PIPE, "/usr/bin/h5dump -d \"$label\" \"$file\" |") or die "$0: Cannot open pipe : $!\n";
	while (<PIPE>) {
		chomp;
		if (/\([0-9]+\)\: +(.+$)/) {
			$value .= $1;
		}
	}

	close(PIPE);

	die "$0: $label not found in $file" if ($value eq "");
	my @a = split/,/, $value;
	trimSpacesVector(\@a);
	return (scalar(@a) == 1) ? $a[0] : @a;
}

sub trimSpacesVector
{
	my ($a) = @_;
	my $n = scalar(@$a);
	for (my $i = 0; $i < $n; ++$i) {
		$a->[$i] =~ s/^ +//;
		$a->[$i] =~ s/ +$//;
	}
}

sub printVector
{
	my ($a) = @_;
	my $n = scalar(@$a);
	for (my $i = 0; $i < $n; ++$i) {
		print "$a->[$i]\n";
	}
}

