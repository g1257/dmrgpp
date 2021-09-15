#!/usr/bin/perl

use strict;
use warnings;
use utf8;

my ($file, $momentumK) = @ARGV;
defined($momentumK) or die "USAGE: $0 filename momentumK\n";

my $label1 = "site <Q0|c'|Q1> time";

my @cdaggerI = loadVector($file, $label1);

my $label2 = "<Q0|c;c'|Q0>";

my @cIcDaggerJ = loadMatrix($file, $label2);

my @numerator = ftVector(\@cdaggerI, $momentumK);

print "$0: numerator= @numerator\n";

my @den = ftMatrix(\@cIcDaggerJ, $momentumK);

print "$0: denominator= @den\n";

die "$0: Denominator isn't real\n" if (abs($den[1]) > 1e-5);

print "$0: Z = ".$numerator[0]/$den[0]." ".$numerator[1]/$den[0]."\n";

sub ftMatrix
{
	my ($a, $momentumK) = @_;
	my $n = scalar(@$a);
	my $center = $n/2;
	my ($sumR, $sumI) = (0, 0);
	for (my $i = 0; $i < $n; ++$i) {
		my $ptr = $a->[$i];
		my $m = scalar(@$ptr);
		($n == $m) or die "$0: Matrix not square\n";
		for (my $j = 0; $j < $m; ++$j) {
			my $arg = $momentumK*($i - $j);
			my $value = ($i < $j) ? $a->[$i]->[$j] : $a->[$j]->[$i];
			$sumR += $value*cos($arg);
			$sumI += $value*sin($arg);
		}
	}

	return ($sumR, $sumI);
}

sub ftVector
{
	my ($a, $momentumK) = @_;
	my $n = scalar(@$a);
	my $center = $n/2;
	my ($sumR, $sumI) = (0, 0);
	for (my $i = 0; $i < $n; ++$i) {
		my $arg = $momentumK*($i - $center);
		my $value = $a->[$i];
		$sumR += $value*cos($arg);
		$sumI += $value*sin($arg);
	}
	
	return ($sumR, $sumI);
}


sub loadVector
{
	my ($file, $label) = @_;
	open(FILE, "<", $file) or die "$0: Cannot open $file : $!\n";
	while (<FILE>) {
		last if (/^\Q$label/);
	}
	
	my @a;
	while (<FILE>) {
		chomp;
		next if (/^#/);
		my @temp = split;
		last if (scalar(@temp) != 3);
		last unless ($temp[0] =~ /^\d+$/);
		$a[$temp[0]] = $temp[1];
	}
	
	close(FILE);
	print STDERR "$0: Found ".scalar(@a)." sites for \"$label\" in $file\n";
	return @a;
}

sub loadMatrix
{
	my ($file, $label) = @_;
	open(FILE, "<", $file) or die "$0: Cannot open $file : $!\n";
	while (<FILE>) {
		last if (/^\Q$label/);
	}
	
	$_ = <FILE>;
	chomp;
	my @sizes = split;
	die "$0: No matrix dimensions found for $label in $file\n" unless (scalar(@sizes) == 2);
	
	my @a;
	my ($rows, $cols) = @sizes;
	for (my $i = 0; $i < $rows; ++$i) {
		$_ = <FILE>;
		chomp;
		next if (/^#/);
		my @temp = split;
		(scalar(@temp) == $cols) or die "$0: row $i not with $cols columns\n";
		$a[$i] = \@temp;
	}
	
	close(FILE);
	print STDERR "$0: Found matrix with ".scalar(@a)." rows for \"$label\" in $file\n";
	return @a;
}


