#!/usr/bin/perl

use strict;
use warnings;
use utf8;

my ($file, $root) = @ARGV;
defined($file) or die "USAGE: $0 filename [outputname]\n";

my @labels = qw/Sigma SiteExcludedG LatticeG/;

defined($root) or $root = getBasename($file);

foreach my $label (@labels) {
	my $buffer = extract($label, $file);
	my $fout = $root.$label;
	writeBuffer($fout, $buffer);
}

sub extract
{
	my ($label, $file) = @_;

	my $flag = 0;
	open(FILE, "<", "$file") or die "$0: Cannot open $file : $!\n";
	while (<FILE>) {
		if (/^$label/) {
			$flag = 1;
			last;
		}
	}

	die "$0: Could not find $label in $file\n" if (!$flag);

	my $size = <FILE>;
	chomp($size);
	print STDERR "$0: Found size=$size for $label in $file\n";

	my $buffer = "";
	for (my $i = 0; $i < $size; ++$i) {
		$_ = <FILE>;
		die "$0: File $file ended before $size elements where read\n" unless ($_);
		chomp;
		my @temp = split;
		die "$0: File $file INTERNAL ERROR for $label\n" if (scalar(@temp) != 2);
		my $omega = $temp[0];
		my $value = $temp[1];
		$value =~ s/\(//;
		$value =~ s/\)//;
		$value =~ s/,/ /;
		$buffer .= "$omega $value\n";
	}

	close(FILE);
	return $buffer;
}

sub writeBuffer
{
	my ($fout, $buffer) = @_;
	open(FOUT, ">", "$fout") or die "$0: Cannot open $fout for writing : $!\n";
	print FOUT "$buffer";
	close(FOUT);
}

