#!/usr/bin/perl

use strict;
use warnings;

package CollectBrakets;

sub main
{
	my ($filename, $foutname) = @_;
	open(FILE, "<", $filename) or die "$0: Cannot open $filename : $!\n";
	if (!open(FOUT, ">", $foutname)) {
		close(FILE);
		die "$0: Cannot write to $foutname : $!\n";
	}

	while (<FILE>) {
		next if (/^#/);
		my $line = $_;
		chomp;
		my @temp = split;
		my $n = scalar(@temp);
		next if ($n != 5);
		next unless isBraket($temp[3]);
		print FOUT $line;
	}

	close(FOUT);
	close(FILE);
}

sub isBraket
{
	my ($x) = @_;
	my ($bra, $ket);
	if ($x =~ /^\<(.+)\|.+\|(.+)\>$/) {
		$bra = $1;
		$ket = $2;
	} else {
		return 0;
	}

	(defined($bra) and defined($ket)) or return 0;

	return (isState($bra) and isState($ket));
}

sub isState
{
	my ($s) = @_;
	return ($s eq "gs" or $s =~ /^P\d+/);
}

1;

