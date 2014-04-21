#!/usr/bin/perl

use warnings;
use strict;

package Utils;

sub getLabel
{
	my ($file,$label)=@_;
	open(FILE,"$file") or die "$0: Cannot open $file: $!\n";
	my $value;
	while(<FILE>) {
		chomp;
		if (/$label(.*$)/) {
			$value=$1;
			last;
		}
	}

	close(FILE);

	defined($value) or die "$0: undefined $label in $file\n";

	return $value;
}

sub myswap
{
	my ($x,$y) = @_;
	my $tmp = $$x;
	$$x = $$y;
	$$y = $tmp;
}

sub setVector
{
	my ($a,$value) = @_;
	my $total = scalar(@$a);
	for (my $i = 0; $i < $total; $i++) {
		$a->[$i] = $value;
	}
}

sub checkRange
{
	my ($val,@list) = @_;

	foreach my $item (@list) {
		return if ($val eq $item);
	}

	die "$0: checkRange: $val does not belong to @list\n";
}

sub loadParams
{
	my ($file) = @_;
	my %params;
	open(FILE,"$file") or die "$0: Cannot open $file\n";
	while(<FILE>) {
		next if (/^#/);
		if (/(^[a-zA-Z]+)=(.*$)/) {
			$params{"$1"}=$2;
		}
	}

	close(FILE);

	return %params;
}

sub scale
{
	my ($file,$f,$c) = @_;

	open(FILE,"$file") or die "$0: Cannot open $file: $!\n";

	my $counter = 0;
	my @saved;
	while(<FILE>) {
		my @temp = split;
		my $n = scalar(@temp);
		if ($n == $c) {
			for (my $i = 1; $i < scalar(@temp); ++$i) {
				$temp[$i] *= $f;
			}
		}

		$saved[$counter++] = \@temp;
	}

	close(FILE);

	open(FOUT,">$file") or die "$0: Cannot open $file: $!\n";

	for (my $i = 0; $i < $counter; ++$i) {
		my $temp = $saved[$i];
		print FOUT "@$temp\n";
	}

	close(FOUT);
}

1;

