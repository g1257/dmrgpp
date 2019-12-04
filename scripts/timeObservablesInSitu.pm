#!/usr/bin/perl

use strict;
use warnings;

package timeObservablesInSitu;

sub getMatrix 
{
	my ($site, $label, $fin, $tau)=@_;
	defined($tau) or return;

	my $counter=0;

	while (<$fin>) {

		last if (/^ALL OPERATORS/);
	}

	my @a;
	while (<$fin>) {
		chomp;
		next unless /\Q$label/;

		if (/^${site} /) {
			my @temp = split;
			(scalar(@temp) >= 5) or next;
			my $value = procValue($temp[1]);
			my $time = int($temp[2]/$tau);
			my $superdensity = procValue($temp[4]);
			my $h = { "value" => $value, "superdensity" => $superdensity};
			$a[$time] = $h;
		}
	}

	return @a;
}

sub load
{
	my ($file) = @_;
	my %h;
	defined($file) or return %h;
	if (!open(FILE, "$file")) {
		print "$0: Could not open $file : $!\n";
	 	return %h;
	}

	$_ = <FILE>;
	defined($_) or return %h;
	chomp;
	if (/#site= (.*$)/) {
		$h{"site"} = $1;
	} else {
		return %h;
	}

	$_ = <FILE>;
	defined($_) or return %h;
	chomp;
	if (/#label=(.*$)/) {
		$h{"label"} = $1;
	} else {
		return %h;
	}

	my $cols = 2;
	my $row = 0;
	my @t;
	my @m;
	while (<FILE>) {
		chomp;
		my @temp = split;
		(scalar(@temp) == 3) or next;
		$t[$row] = $temp[0];
		$m[0 + $row*$cols] = $temp[1];
		$m[1 + $row*$cols] = $temp[2];
		++$row;
	}

	close(FILE);

	$h{"times"} = \@t;
	my $rows = $row;
	my @m2;
	$m2[0] = $rows;
	$m2[1] = $cols;
	for (my $i = 0; $i < $rows; ++$i) {
		for (my $j = 0; $j < $cols; ++$j) {
			$m2[2 + $i + $j*$rows] = $m[$j + $i*$cols];
		}
	}

	$h{"data"} = \@m2;
	return %h;
}

sub procValue
{
	my ($t)=@_;
	$_=$t;
	s/\(//;
	s/,.*$//;
	return $_;
}

1;
