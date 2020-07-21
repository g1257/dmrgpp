#!/usr/bin/perl

use strict;
use warnings;
use utf8;
use Math::Trig;

my $pi = Math::Trig::pi;
my ($file, $label) = @ARGV;
defined($label) or die "USAGE: $0 filename label\n";

my @m = loadMatrix($file, $label);

my @ft = fourierTransform(\@m);

printArray(\@ft);

sub printArray
{
	my ($array) = @_;
	my $n = scalar(@$array);
	for (my $m = 0; $m < $n; ++$m) {
		my $q = 2*pi*$m/$n;
		print "$q ".$array->[$m]."\n";
	}
}

sub fourierTransform
{
	my ($array) = @_;
	my $n = scalar(@$array);
	my @sq;
	my $c = int($n/2);
	for (my $m = 0; $m < $n; ++$m) {
		
		my $sum = 0;
		
		for (my $i = 0; $i < $n; ++$i) {
			my $ptr1 = $array->[$i];
			my $ptr2 = $array->[$c];
			die "$0: undefined for $i \n" if (!defined($ptr1));
			die "$0: undefined for $c \n" if (!defined($ptr2));

			
			my $value = ($i < $c) ? $array->[$i]->[$c] : $array->[$c]->[$i];
			
			$sum += cos(2*pi*($i-$c)*$m/$n)*$value;
		}
		
		$sq[$m] = $sum;
	}
	
	return @sq;
}

sub loadMatrix
{
	my ($file, $label) = @_;
	open(FILE, "<", "$file") or die "$0: Cannot open $file : $!\n";
	my @m;
	my $found = 0;
	while (<FILE>) {
		next if (/cmdline:/i);
		if (/^$label/) {
			$found = 1;
			last;
		}
	}
	
	if (!$found) {
		close(FILE);
		die "$0: Cannot find $label in $file\n" 
	}
	
	$_ = <FILE>;
	chomp;
	my @temp = split;
	if (scalar(@temp) != 2) {
		close(FILE);
		die "$0: Expected 2 sizes found @temp instead\n";
	}
	
	my $rows = $temp[0];
	my $cols = $temp[1];
	
	while (<FILE>) {
		chomp;
		my @temp = split;
		last if (scalar(@temp) != $cols);
		push @m, \@temp;
	}
	
	close(FILE);
	my $r = scalar(@m);
	if ($r != $rows) {
		die "$0: Expected $rows rows but found $r instead\n";
	}
	 
	
	return @m;
}

