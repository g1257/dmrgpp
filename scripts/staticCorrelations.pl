#!/usr/bin/perl

use strict;
use warnings;
use utf8;

my ($file, $label) = @ARGV;

my @m = readMatrix($file, $label);

@m = symmetrize(@m);

my $rows = scalar(@m);

for (my $i = 0; $i < $rows; ++$i) {
	my $ptr = $m[$i];
	my $cols = scalar(@$ptr);
	for (my $j = 0; $j < $cols; ++$j) {
		print $ptr->[$j]." ";
	}
	
	print "\n";
}

sub symmetrize
{
	my @m = @_;
	my @m2;
	my $rows = scalar(@m);
	for (my $i = 0; $i < $rows; ++$i) {
		my $ptr = $m[$i];
		my $cols = scalar(@$ptr);
		my @tmp = @$ptr;
		for (my $j = 0; $j < $i; ++$j) {
				@tmp[$j] = $m[$j]->[$i];
		}
	
		$m2[$i] = \@tmp;
}

		return @m2;
}


sub readMatrix
{
	my ($file, $label) = @_;
	open(FILE, "<", $file) or die "$0: Cannot open $file : $!\n";
	while (<FILE>) {
		next if (/CmdLine/);
		last if (/$label/);
	}

	$_ = <FILE>;
	chomp;
	my @temp = split(/\s/, $_);
	die "$0: Rows Cols not found, line $. ** $_ **\n" if (scalar(@temp) != 2);
	my $rows = $temp[0];
	my $cols = $temp[1];

	my $row = 0;
	my @m;
	while (<FILE>) {
		chomp;
		my @tmp = split;
		die "$0: Wrong number of cols, line $., $cols not equal "
		.scalar(@tmp)."\n" if ($cols != scalar(@tmp));

		$m[$row++] = \@tmp;
		last if ($row == $rows);
	}

	return @m;
}
