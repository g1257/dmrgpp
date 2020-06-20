#!/usr/bin/perl

use strict;
use warnings;
use utf8;
use Math::Trig;

my ($file, $wbegin, $wtotal, $wstep, $mode) = @ARGV;
defined($wstep) or die "USAGE: $0 filename wbegin wtotal wstep [mode]\n";
defined($mode) or $mode = -1;
my @data;

loadData(\@data, $file);
my $n = scalar(@data);
print STDERR "$0: Found $n sites\n";
printLoadedData(\@data, $mode) if ($mode >= 0);

my @omegas = fillOmegas($wbegin, $wtotal, $wstep);
my @dataOmega = ftTime(\@data, \@omegas);
printSpaceOmega(\@dataOmega, \@omegas) if ($mode == -2);

my @dataOmegaK = ftSpace(\@dataOmega, scalar(@omegas));

printData(\@dataOmegaK, \@omegas) if ($mode == -1);

sub loadData
{
	my ($data, $file) = @_;
	open(FILE, "<", "$file") or die "$0: Cannot open $file : $!\n";
	while (<FILE>) {
		next if (/^#/);
		my @temp = split;
		next if (scalar(@temp) != 4);
		my $time = $temp[0];
		my $site = $temp[1];
		my $value = $temp[2];
		$value = 0 if ($value == -100);
		my $h = {"time" => $time, "value" => $value};
		
		if (defined($data->[$site])) {
			my $a = $data->[$site];
			push @$a, $h;
		} else {
			my @a = ($h);
			$data->[$site] = \@a;
		}
	}

	close(FILE);
}

sub printLoadedData
{
	my ($data, $mode) = @_;
	my $start = 0;
	my $end = scalar(@$data);
	if ($mode >= 0 && $mode < $end) {
		$start = $mode;
		$end = $start + 1;
	}

	my $ptr = $data->[1];
	my $nptr = scalar(@$ptr);
	my %hash;
	for (my $i = 0; $i < $nptr; ++$i) {
		my $h = $ptr->[$i];
		my $time = $h->{"time"};
		my $value = $h->{"value"};
		$hash{$time} = $value;
	}

	for my $time (sort keys %hash) {
		print "$time ";
		for (my $i = $start; $i < $end; ++$i) {
			my $ptr = $data->[$i];
			my $nptr = scalar(@$ptr);
			for (my $j = 0; $j < $nptr; ++$j) {
				my $h = $ptr->[$j];
				my $t = $h->{"time"};
				next unless ($time == $t);
				my $value = $h->{"value"};
				print " $value";
			}
		}

		print "\n";
	}
}

sub fillOmegas
{
	my ($wbegin, $wtotal, $wstep) = @_;
	my @omegas;
	for (my $i = 0; $i < $wtotal; ++$i) {
		my $omega = $wbegin + $i*$wstep;
		push @omegas, $omega;
	}

	return @omegas;
}

sub ftTime
{
	my ($data, $omegas) = @_;
	my @dataO;
	for (my $i = 0; $i < $n; ++$i) {
		my $src = $data->[$i];
		my @dest = ftTimeOneSite($src, $omegas);
		$dataO[$i] = \@dest;
	}

	return @dataO;
}

sub printSpaceOmega
{
	my ($dataOmega, $omegas) = @_;
	my $n = scalar(@omegas);
	my $sites = scalar(@dataOmega);
	my $center = int($sites/2);
	for (my $i = 0; $i < $n; ++$i) {
		my $omega = $omegas->[$i];
		my $value = $dataOmega->[$center]->[$i];
		print "$omega $value\n";
	}
}


sub ftTimeOneSite
{
	my ($src, $omegas) = @_;
	# a(omega) = \sum_times sin(omega*t) * src[t]
	my @dest;
	my $n = scalar(@$omegas);
	my $m = scalar(@$src);
	for (my $i = 0; $i < $n; ++$i) {
		my $omega = $omegas->[$i];
		my $sum = 0;
		for (my $j = 0; $j < $m; ++$j) {
			my $ptr = $src->[$j];
			my $time = $ptr->{"time"};
			my $value = $ptr->{"value"};
			$sum += sin($omega*$time)*$value;
		}

		$dest[$i] = $sum;
	}

	return @dest;
}

sub ftSpace
{
	my ($src, $omegas) = @_;
	
	my @dest;
	my $ks = scalar(@$src);
	my $sites = $ks;
	for (my $k = 0; $k < $ks; ++$k) {
		my @dummy;
		$dest[$k] = \@dummy;
	}
 
	for (my $i = 0; $i < $omegas; ++$i) {
		my @values;
		for (my $j = 0; $j < $sites; ++$j) {
			$values[$j] = $src->[$j]->[$i];
		}

		my @tmp = ftSpaceOneOmega(\@values);
	
		for (my $k = 0; $k < $ks; ++$k) {
			$dest[$k]->[$i] = $tmp[$k];
		}
	}

	return @dest;
}

sub printData
{
	my ($vals, $omegas) = @_;
	my $pi = Math::Trig::pi;
	my $n = scalar(@$omegas);
	my $sites = scalar(@$vals);
	for (my $i = 0; $i < $n; ++$i) {
		my $omega = $omegas->[$i];
		for (my $k = 0; $k < $sites; ++$k) {
			my $kactual = 2*$k*$pi/$sites;
			my $value = $vals->[$k]->[$i];
			$value = 0 if (fabs($value) < 1e-4);
			$value = int($value*1000)/1000;
			print "$kactual $omega $value\n";
		}

		print "\n";
	}
}

sub fabs
{
	my ($x) = @_;
	return ($x >= 0) ? $x : -$x;
}

sub ftSpaceOneOmega
{
	my ($array) = @_;
	my $n = scalar(@$array);
	my $center = int($n/2);
	my @dest;
	my $pi = Math::Trig::pi;
	for (my $k = 0; $k < $n; ++$k) {
		my $kactual = 2*$k*$pi/$n;
		my $sum = 0;
		for (my $i = 0; $i < $n; ++$i) {
			$sum += $array->[$i]*cos($kactual*($i-$center));
		}

		$dest[$k] = $sum;
	}

	return @dest;
}


