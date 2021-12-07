#!/usr/bin/perl

use strict;
use warnings;
use utf8;

my @models = @ARGV;

my $file = "Engine/ModelSelector.h0";
my $foutname = "Engine/ModelSelector.h";

if (scalar(@models) == 0) {
	system("cp $file $foutname");
	print STDERR "Restored $foutname\n";
	exit(0);
}

my $finput;
open($finput, "<", $file) or die "$0: Cannot open $file : $!\n";
my $fout;
open($fout, ">", $foutname) or die "$0: Cannot write to $foutname : $!\n";

readAndWriteUntil($fout, "// start headers", $finput);

readAndWriteLinesThatMatch($fout, \@models, $finput);

readAndWriteUntil($fout, "// start models here", $finput);

readAndWriteLinesThatMatch($fout, \@models, $finput);

readAndWriteUntil($fout, "// named models start", $finput);

readAndWriteLinesThatMatch($fout, \@models, $finput, 1);

readAndWriteUntil($fout, "", $finput);
close($fout);
close($finput);

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
	my ($fout, $a, $fh, $isNames) = @_;
	defined($isNames) or $isNames = 0;
	my $count = 0;
	while (<$fh>) {
		my $line = $_;
		last if (/\/\/ end models/);

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

