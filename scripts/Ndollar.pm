#!/usr/bin/perl

use strict;
use warnings;
use utf8;

package Ndollar;

sub main
{
	my ($file) = @_;
	defined($file) or die "USAGE: $0 file\n";

	my $label = "site <gs|n|gs> time";

	open(FILE, "<", $file) or die "$0: Cannot open $file : $!\n";
	while (<FILE>) {
		last if (/\Q$label/);
	}

	my @values;

	while (<FILE>) {
		chomp;
		my @temp = split;
		(scalar(@temp) == 3) or next;
		my $site = $temp[0];
		($site =~ /^\d+$/) or next;
		$values[$site] = $temp[1];
	}

	close(FILE);

	my $n = scalar(@values);

	print STDERR "$0: Read $n values\n";

	die "$0: Too big $n >= 100\n" if ($n >= 100);

	for (my $i = 0; $i < $n; ++$i) {
		my $v0 = -$values[$i];
		my $v1 = 1.0 + $v0;
		my $v2 = 2.0 + $v0;
		my $fout = "n$i.txt";
		open(FOUT, ">", "$fout") or die "$0: Cannot write to $fout : $!\n"; 
		print FOUT<<EOF;
TSPOperator=raw
RAW_MATRIX
4 4
$v0 0 0 0
0 $v1 0 0
0 0 $v1 0
0 0 0 $v2
FERMIONSIGN=1
JMVALUES 2 0 0
AngularFactor=1

EOF
		close(FOUT);
	}

	print STDERR "$0: Written n*.txt $n files\n";
}

1;


