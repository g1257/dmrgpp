#!/usr/bin/perl

use strict;
use warnings;
use utf8;

package MinCompile;

sub main
{
	my ($what, $array) = @_;

	if ($what ne "models" and $what ne "targets") {
		die "$0: Expected models or targets not $what\n";
	}

	my $capital = $what;
	$capital =~ s/s$//;
	$capital =~ s/(^.{1})/\U$1/;
	my $file = "Engine/".$capital."Selector.h0";
	my $foutname = "Engine/".$capital."Selector.h";

	if (scalar(@$array) == 0) {
		system("cp $file $foutname");
		print STDERR "Restored $foutname\n";
		return 0;
	}

	my $finput;
	open($finput, "<", $file) or die "$0: Cannot open $file : $!\n";
	my $fout;
	open($fout, ">", $foutname) or die "$0: Cannot write to $foutname : $!\n";

	readAndWriteUntil($fout, "// start headers", $finput);

	readAndWriteLinesThatMatch($fout, $array, $finput, $what);

	readAndWriteUntil($fout, "// start $what here", $finput);

	readAndWriteLinesThatMatch($fout, $array, $finput, $what);

	readAndWriteUntil($fout, "// named $what start", $finput);

	my $stripElse = ($what eq "models");
	readAndWriteLinesThatMatch($fout, $array, $finput, $what, $stripElse);

	readAndWriteUntil($fout, "", $finput);
	close($fout);
	close($finput);
	print STDERR "$0: File $foutname has been written.\n";
}

sub readAndWriteUntil
{
	my ($fout, $text, $fh) = @_;
	my $count = 0;
	while (<$fh>) {
		last if (/$text/ and $text ne "");
		print {$fout} $_;
		++$count;
	}
}

sub readAndWriteLinesThatMatch
{
	my ($fout, $a, $fh, $what, $isNames) = @_;
	defined($isNames) or $isNames = 0;
	my $count = 0;
	while (<$fh>) {
		my $line = $_;
		last if (/\/\/ end $what/);

		if (isInArray($line, $a)) {
			s/\} else // if ($count++ == 0 and $isNames);
			print {$fout} $_;
		}
	}
}

sub isInArray
{
	my ($line, $a) = @_;
	my $n = scalar(@$a);
	for (my $i = 0; $i < $n; ++$i) {
		my $val = $a->[$i];
		return 1 if ($line =~ /$val/);
	}

	return 0;
}

1;
