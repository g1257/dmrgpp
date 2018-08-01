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
	die "$0: File $file not readable\n" unless (-r "$file");
	die "$0: Label $label is invalid\n" unless isValidLabel($label);

	my $value = "";
	open (PIPE, "h5dump -d \"$label\" \"$file\" |") or die "$0: Cannot open pipe : $!\n";
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

sub isValidLabel
{
	my ($label) = @_;
	my @a = split/\//, $label;
	my $n = scalar(@a);
	for (my $i = 0; $i < $n; ++$i) {
		return 0 unless isValidSubLabel($a[$i]);
	}

	return 1;
}

sub isValidSubLabel
{
	my ($slabel) = @_;
	return ($slabel eq "" or $slabel =~ /^[a-zA-Z0-9]+$/) ? 1 : 0;
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

