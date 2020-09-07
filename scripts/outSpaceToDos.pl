#!/usr/bin/perl

use strict;
use warnings;
use utf8;

my ($file, $site) = @ARGV;
defined($site) or die "USAGE: $0 file site\n";

my %data = readData($file);

plotData(\%data, $site);

sub readData
{
	my ($file) = @_;
	open(FILE, "<", $file) or die "$0: Cannot open file : $!\n";

	my $n;
	my %data;
	while (<FILE>) {
		chomp;
		my @temp = split; # omega, n
		my $x = scalar(@temp);
		die "$0: FATAL: $file : $. : Expecting 2 items, got $n instead\n" if ($x != 2);
		my $omega = $temp[0];
		if (!defined($n)) {
			$n = $temp[1];
		} else {
			die "$0: FATAL: $file : $., expecting second number to be $n, got ".$temp[1]."\n" if ($n != $temp[1]);
		}

		my @values;
		for (my $i = 0; $i < $n; ++$i) {
			$_ = <FILE>;
			chomp;
			my @temp = split; # site, real, imag
			my $x = scalar(@temp);
			die "$0: FATAL: $file : $. : Expecting 3 items, got $x instead\n" if ($x != 3);
			push @values, \@temp;
		}

		$data{"$omega"} = \@values;
	}

	close(FILE);
	return %data;
}

sub plotData
{
	my ($h, $site) = @_;
	my $hasParens = 0;
	foreach my $omega (sort {$a <=> $b} keys %$h) {
		my $ptr = $h->{$omega};
		my $n = scalar(@$ptr);
		die "FATAL: $site >= $n\n" if ($site >= $n);
		my $ptr2 = $ptr->[$site];
		die "$0: Internal error, expecting 3 elements\n" if (scalar(@$ptr2) != 3);
		my $x = $ptr2->[0];
		die "$0: Internal error, expecting 1st element to be $site, got $x\n" if ($x != $site);
		my $value1 = $ptr2->[1];
		my $value2 = $ptr2->[2];
		$hasParens = checkParens($value1, $value2);
		if ($hasParens) {
			($value1, $value2) = reIm($value2);
		}

		print "$omega $value1 $value2\n";
	}
}

sub checkParens
{
	my ($val1, $val2) = @_;
	my $flag = hasParens($val1);
	die "$0: INTERNAL error (parens)\n" if ($flag && !hasParens($val2));
	return $flag;
}

sub hasParens
{
	my ($val) = @_;
	return ($val =~ /\(/);
}

sub reIm
{
	my ($val) = @_;
	my @temp = split/,/, $val;
	scalar(@temp) == 2 or die "$0: Not a complex value $val\n";
	my $re = $temp[0];
	$re =~ s/\(//;
	my $im = $temp[1];
	$im =~ s/\)//;
	return ($re, $im);
}

