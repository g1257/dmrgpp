#!/usr/bin/perl

use strict;
use utf8;
use warnings;

my ($file) = @ARGV;
defined($file) or die "USAGE: $0 filename\n";

my @entropy = load($file);

printEntropy(\@entropy);

sub printEntropy
{
	my ($entropy) = @_;
	my $n = scalar(@$entropy);
	print STDERR "$0: $n values in array\n";
	my $nOver2 = int($n/2);
	for (my $i = 0; $i < $n; ++$i) {
		my $value = $entropy->[$i];
		defined($value) or next;
		$value = $entropy->[$n - $i] if ($i >= $nOver2);
		print "$i $value\n";
	}
}

sub load
{
	my ($file) = @_;
	open(FILE, "<", $file) or die "$0: Cannot open $file : $!\n";
	my $split;
	my $flag = 0;
	my @entropy;
	while (<FILE>) {
		if (/sites=([^\+]+)\+/) {
			$split = $1;
			$flag = 1;
		}
		
		if (/EntropyVonNeumann= ([^;]+);/) {
			my $entropy = $1;
			#last if ($flag != 1);
			$flag = 0;
			$entropy[$split] = $entropy;
		}	
	}
	
	close(FILE);
	
	return @entropy;
}
	
