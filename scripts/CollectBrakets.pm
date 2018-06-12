#!/usr/bin/perl

use strict;
use warnings;

package CollectBrakets;

sub main
{
	my ($filename, $foutname) = @_;
	open(FILE, "<", $filename) or die "$0: Cannot open $filename : $!\n";
	if (!open(FOUT, ">", $foutname)) {
		close(FILE);
		die "$0: Cannot write to $foutname : $!\n";
	}

	while (<FILE>) {
		next if (/^#/);
		my $line = $_;
		chomp;
		my @temp = split;
		my $n = scalar(@temp);
		next if ($n != 5);
		next unless isBraket($temp[3]);
		print FOUT almostZeroToZero($line);
	}

	close(FOUT);
	close(FILE);
}

sub almostZeroToZero
{
	my ($x) = @_;
	$_ = $x;
	my $hasEol = 0;
	if (substr($x, -1) eq "\n") {
		$hasEol = 1;
	} else {
		$_ = $x;
	}

	my @temp = split;
	my $ret = "";
	my $n = scalar(@temp);
	for (my $i = 0; $i < $n; ++$i) {
		$ret .= almost0To0($temp[$i])." ";
	}

	return ($hasEol) ? $ret."\n" : $ret;
}

sub almost0To0
{
	my ($x) = @_;
	return $x if (!isFloat($x));
	return (abs($x) < 1e-6) ? "0" : $x;
}

sub isFloat
{
	my ($x) = @_;
	$_ = $x;
	return (/^[+-]?(?=\.?\d)\d*\.?\d*(?:e[+-]?\d+)?\z/i);
}

sub isBraket
{
	my ($x) = @_;
	my ($bra, $ket);
	if ($x =~ /^\<(.+)\|.+\|(.+)\>$/) {
		$bra = $1;
		$ket = $2;
	} else {
		return 0;
	}

	(defined($bra) and defined($ket)) or return 0;

	return (isState($bra) and isState($ket));
}

sub isState
{
	my ($s) = @_;
	return ($s eq "gs" or $s =~ /^P\d+/);
}

1;

