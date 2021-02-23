#!/usr/bin/perl

use strict;
use warnings;
use utf8;

my ($file, $n) = @ARGV;
defined($file) or die "USAGE: $0 filename [set number]\n";
defined($n) or $n = -1;

my @data = loadData($file, $n);

printData(\@data);

sub printData
{
	my ($a) = @_;
	my $n = scalar(@$a);
	for (my $i = 0; $i < $n; ++$i) {
		print "$i ".$a->[$i]."\n";
	}
}

sub loadData
{
	my ($file, $n) = @_;
	open(FILE, "<", "$file") or die "$0: Cannot open $file : $!\n";
	my @temp;

	my $counter = 0;
	while (<FILE>) {
		next unless (/RitzEigenvalues: index\=\d+\|/);
		s/RitzEigenvalues: index=\d+\|//;
		chomp;
		@temp = split;
		last if ($counter++ == $n);
	}

	close(FILE);
	my $m = scalar(@temp);
	die "$0: No /RitzEigenvalues found in $file\n" if ($m == 0);
	return @temp;
}


