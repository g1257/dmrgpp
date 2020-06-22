#!/usr/bin/perl

use strict;
use warnings;
use utf8;
use Math::Trig;

my $pi = Math::Trig::pi;
my ($file, $wbegin, $wtotal, $wstep, $mode) = @ARGV;
defined($wstep) or die "USAGE: $0 filename wbegin wtotal wstep [mode]\n";
defined($mode) or $mode = -1;
my @data;

my $tmax = loadData(\@data, $file);
my $n = scalar(@data);
print STDERR "$0: Found $n sites\n";
printLoadedData(\@data, $mode) if ($mode >= 0);

my @omegas = fillOmegas($wbegin, $wtotal, $wstep);
my @dataOmega = ftTime(\@data, \@omegas, \&dampFunction);
printSpaceOmega(\@dataOmega, \@omegas) if ($mode == -2);

my @dataOmegaK = ftSpace(\@dataOmega, scalar(@omegas));

printData(\@dataOmegaK, \@omegas) if ($mode == -1);

sub dampFunction
{
	my ($time) = @_;
	return 0.5*(1.0+cos(($time)*pi/$tmax));
}

sub loadData
{
	my ($data, $file) = @_;
	my $tmax = 0;
	open(FILE, "<", "$file") or die "$0: Cannot open $file : $!\n";
	while (<FILE>) {
		next if (/^#/);
		my @temp = split;
		next if (scalar(@temp) != 4);
		my $time = $temp[0];
		my $site = $temp[1];
		my $value = $temp[2];
		$value = "(0,0)" if ($value eq "-100");
		my $h = {"time" => $time, "value" => $value};

		$tmax = $time if ($time > $tmax);

		if (defined($data->[$site])) {
			my $a = $data->[$site];
			push @$a, $h;
		} else {
			my @a = ($h);
			$data->[$site] = \@a;
		}
	}

	close(FILE);
	return $tmax;
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
	my ($data, $omegas, $dampF) = @_;
	my @dataO;
	for (my $i = 0; $i < $n; ++$i) {
		my $src = $data->[$i];
		my @dest = ftTimeOneSite($src, $omegas, $dampF);
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
		my ($re, $im) = realImag($dataOmega->[$center]->[$i]);
		print "$omega $re $im\n";
	}
}


sub ftTimeOneSite
{
	my ($src, $omegas, $dampF) = @_;
	# a(omega) = \sum_times sin(omega*t) * src[t]
	my @dest;
	my $n = scalar(@$omegas);
	my $m = scalar(@$src);
	for (my $i = 0; $i < $n; ++$i) {
		my $omega = $omegas->[$i];
		my ($sumr, $sumi) = (0, 0);
		for (my $j = 0; $j < $m; ++$j) {
			my $ptr = $src->[$j];
			my $time = $ptr->{"time"};
			my $value = $ptr->{"value"};
			my ($c, $s) = (cos($omega*$time), sin($omega*$time));
			my ($re, $im) = realImag($value);
			my $damp = $dampF->($time);
			$sumr += ($c*$re - $s*$im)*$damp;
			$sumi += ($c*$im + $s*$re)*$damp;
		}

		$dest[$i] = "(".$sumr.",".$sumi.")";
	}

	return @dest;
}

sub realImag
{
	my ($x) = @_;
	my @temp = split/,/, $x;
	my $n = scalar(@temp);
	$n == 2 or die "$0: Complex number $x\n";
	my $re = $temp[0];
	my $im = $temp[1];
	$re =~ s/\(//;
	$im =~ s/\)//;
	return ($re, $im);
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
	my $n = scalar(@$omegas);
	my $sites = scalar(@$vals);
	for (my $i = 0; $i < $n; ++$i) {
		my $omega = $omegas->[$i];
		for (my $k = 0; $k < $sites; ++$k) {
			my $kactual = 2*$k*$pi/$sites;
			my $value = $vals->[$k]->[$i];
			my ($re, $im) = realImag($value);
			#$value = 0 if (fabs($value) < 1e-4);
			#$value = int($value*1000)/1000;
			print "$kactual $omega $im\n";
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
		my ($sumr, $sumi) = (0, 0);
		for (my $i = 0; $i < $n; ++$i) {
			my ($re, $im) = realImag($array->[$i]);
			my $arg = $kactual*($i-$center);
			my ($c, $s) = (cos($arg), sin($arg));
			$sumr += $c*$re - $s*$im;
			$sumi += $c*$im + $s*$re;
			#$sum += $array->[$i]*cos($kactual*($i-$center));
		}

		$dest[$k] = "(".$sumr.",".$sumi.")";
	}

	return @dest;
}

