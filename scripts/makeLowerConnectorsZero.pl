#!/usr/bin/perl

use strict;
use warnings;
use utf8;

my ($file) = @ARGV;
defined($file) or die "USAGE: $0 filename\n";

my @matrix3 = readMatrix($file);

for (my $i = 0; $i < scalar(@matrix3); ++$i) {
	my @matrix2 = setLowerPartZero($matrix3[$i]);
	printMatrix(\@matrix2);
}

sub printMatrix
{
	my ($matrix) = @_;
	my $rows = scalar(@$matrix);

	 for (my $i = 0; $i < $rows; ++$i) {
		  my $ptr = $matrix->[$i];
		  my $cols = scalar(@$ptr);
		  for (my $j = 0; $j < $cols; ++$j) {
			  my $val = $ptr->[$j];
			  print "$val ";
		  }

		  print "\n";
	  }

	  print "\n\n\n";
}

sub setLowerPartZero
{
	my ($matrix) = @_;
	my $rows = scalar(@$matrix);
	my @matrix2;
	my $cols;
	for (my $i = 0; $i < $rows; ++$i) {
		my $ptr = $matrix->[$i];
		my $m = scalar(@$ptr);
		if (defined($cols) and $cols != $m) {
			die "$0: setLowerPartZero: cols of different sizes $m $cols\n";
		} else {
			$cols = $m;
		}

		my @temp;
		for (my $j = 0; $j < $cols; ++$j) {
			$temp[$j] = ($i > $j) ? 0 : $ptr->[$j];
		}

		push @matrix2, \@temp;
	}

	return @matrix2;
}

sub readMatrix
{
	my ($file) = @_;
	my $label = "Connectors";
	my @matrix;
	my @matrix3;
	my ($rows, $cols);
	my $flag = 0;
	open(FILE, "<", $file) or die "$0: Cannot open $file : $!\n";
	while (<FILE>) {
		if (/${label}[ \t]+(\d+)[ \t]+(\d+)/) {
			$rows = $1;
			$cols = $2;
			$flag = 1;
			next;
		}

		if ($flag) {
			my @temp = split;
			my $m = scalar(@temp);
			next unless $m > 0;
			if ($m != $cols) {
				die "$0: readMatrix: cols of different sizes $m $cols\n";
			}

			push @matrix, \@temp;
			if (scalar(@matrix) == $rows) {
				my @matrix2 = @matrix;
				push @matrix3, \@matrix2;
				@matrix = ();
				$flag = 0;
				print STDERR "$0: ".scalar(@matrix3)." matrices have been read\n";
				next;
			}

		}
	}

	close(FILE);
	if (!defined($rows) or !defined($cols)) {
		die "$0: Matrix $label not found in $file\n";
	}

	print STDERR "$0: Found matrix with $rows times $cols\n";
	return @matrix3;
}

