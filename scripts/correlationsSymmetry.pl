#!/usr/bin/perl

use strict;
use warnings;

my ($file, $label) = @ARGV;

defined($label) or die "USAGE: $0 file label\n";

my @matrix = readMatrix($file, $label);
my $n = $matrix[0];
my @h = createPairs($n);

#printPairs1(\@h);

examinePairs1(\@matrix, \@h);

sub createPairs
{
	my ($n) = @_;
	my @h;
	my @a;
	my @b;
	my $counter = 0;
	($n & 1) and die "$0: Number of sites must be even\n";
	my $nOver2 = int($n/2);
	for (my $i = 0; $i < $nOver2; ++$i) {
		for (my $j = $i + 1; $j < $n; ++$j) {
			my @a = ($i, $j);
			my @b = findSymmetric($i, $j, $n);
			my @temp = (\@a, \@b);
			$h[$counter++] = \@temp;
		}
	}

	return @h;
}

sub findSymmetric
{
	my ($ind, $jnd, $l) = @_;
	($ind < $jnd) or die "$0: Error at findSymmetric\n";
	($l - 1 >= $ind) or  die "$0: Error at findSymmetric\n";
	my $i = $l - 1 - $ind;
	my $j = $l - 1 - $jnd;
	return ($j, $i);
}

sub printPairs1
{
	my ($h) = @_;
	my $p = scalar(@$h);
	for (my $i = 0; $i < $p; ++$i) {
		my $ptr = $h->[$i];
		(scalar(@$ptr) == 2) or die "$0: bad pairs\n";
		my $a = $ptr->[0];
		my $b = $ptr->[1];
		printPairs2($a, $b);
	}
}

sub printPairs2
{
	my ($a, $b) = @_;
	my @a = @$a;
	my @b = @$b;
	(scalar(@a) == 2) or die "$0: bad pair\n";
	(scalar(@b) == 2) or die "$0: bad pair\n";
	print "($a[0], $a[1]) ==   ($b[0], $b[1])\n";
}

sub readMatrix
{
	my ($file, $label) = @_;
	open(FILE, "$file") or die "$0: Cannot open file $file : $!\n";
	while (<FILE>) {
		last if (/^\Q$label/);
	}

	$_ = <FILE>;
	defined($_) or die "$0: Nothing found after $label in $file\n";
	chomp;
	my @temp = split;
	(scalar(@temp) == 2) or die "$0: Error reading $label in $file\n";
	$matrix[0] = $temp[0];
	$matrix[1] = $temp[1];
	my $rows = $matrix[0];
	my $cols = $matrix[1];
	for (my $i = 0; $i < $rows; ++$i) {
		$_ = <FILE>;
		defined($_) or die "$0: Error reading matrix $label in $file\n";
		chomp;
		my @temp = split;
		(scalar(@temp) == $cols) or die "$0: Error in row $i, $file\n";
		my @t = saveLine(\@temp);
		(scalar(@t) == $cols) or die "$0: Error in row $i, in $file\n";
		for (my $j = 0; $j < $cols; ++$j) {
			$matrix[2 + $j*$rows + $i] = $t[$j];
		}
	}

	close(FILE);
	return @matrix;
}

sub saveLine
{
	my ($a) = @_;
	my $cols = scalar(@$a);
	my @t;
	for (my $i = 0; $i < $cols; ++$i) {
		my $x = $a->[$i];
		$x =~ s/[\(\)]//;
		$x =~ s/,.*$//;
		$t[$i] = $x;
	}

	return @t;
}

sub examinePairs1
{
	my ($m, $h) = @_;
	my $p = scalar(@$h);
	for (my $i = 0; $i < $p; ++$i) {
		my $ptr = $h->[$i];
		(scalar(@$ptr) == 2) or die "$0: bad pairs\n";
		my $a = $ptr->[0];
		my $b = $ptr->[1];
		examinePairs2($m, $a, $b);
	}
}

sub examinePairs2
{
	my ($m, $a, $b) = @_;
	my @a = @$a;
	my @b = @$b;
	(scalar(@a) == 2) or die "$0: bad pair\n";
	(scalar(@b) == 2) or die "$0: bad pair\n";
	my $i = $a[0];
	my $j = $a[1];
	return if ($i == $b[0] and $j == $b[1]);

	my $x = printMatrixElement($m, $i, $j);
	print " --- ";
	$i = $b[0];
	$j = $b[1];
	my $y = printMatrixElement($m, $i, $j);
	my $diff = abs($x - $y);
	print " -- diff = $diff\n";
}		

sub printMatrixElement
{
	my ($m, $i, $j) = @_;
	my $rows = $m->[0];
	my $x = $m->[2 + $j*$rows + $i];
	print "m($i, $j)=".$x;
	return $x;
}

