#!/usr/bin/perl

use strict;
use warnings;
use utf8;
use Fourier;

my ($filename) = @ARGV;
defined($filename) or die "USAGE: $0 filename\n";

my $labeln = "site <gs|n|gs>";
my @n = loadVector($filename, $labeln);
my $averageN = computeAverage(\@n);

my $labelNn = "<gs|n;n|gs>";
my @Nn = Fourier::loadMatrix($filename, $labelNn);

my @nq = Fourier::fourierTransform(\@Nn);

correctNqZero(\@nq, $averageN);

Fourier::printArray(\@nq);

# Here we n_{q=0}Corrected = n_{q=0} -L*n^2, 
sub correctNqZero
{
	my ($nq, $average) = @_;
	my $val = $nq->[0]->[0];
	if (!defined($val)) {
		print STDERR "$0: Could not correct nq[0], it does not exist\n";
		return;
	}
	
	my $nsites = scalar(@$nq);
	my $correction = $nsites*$average*$average;
	$nq->[0]->[0] = $val - $correction;
}

sub loadVector
{
	my ($file, $label) = @_;
	open(FILE, "<", "$file") or die "$0: Cannot open $file : $!\n";
	my $found = 0;
	while (<FILE>) {
		next if (/cmdline:/i);
		if (/^\Q$label/) {
			$found = 1;
			last;
		}
	}

	if (!$found) {
		close(FILE);
		die "$0: Cannot find $label in $file\n"
	}
	
	my @v;
	
	while (<FILE>) {
		my @temp = split;
		last unless (scalar(@temp) == 3);
		$v[$temp[0]] = $temp[1];
	}
	
	close(FILE);
	print STDERR "$0: Read vector $label from $file with ".scalar(@v)." entries.\n";
	return @v;
}

sub computeAverage
{
	my ($v) = @_;
	my $n = scalar(@$v);
	my $sum  = 0;
	for (my $i = 0; $i < $n; ++$i) {
		$sum += $v->[$i];
	}
	
	return $sum/$n;
}


