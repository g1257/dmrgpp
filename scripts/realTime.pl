#!/usr/bin/perl
#
use strict;
use warnings;

my ($input) = @ARGV;
defined ($input) or die "USAGE: $0 input.inp\n";

my $root = $input;

if ($root =~ /\.inp$/) {
	$root =~ s/\.inp$//;
} else {
	$input .= ".inp";
}

my $cout = "runFor$root.cout";
my $t = getRunTime($cout);

print "WallTimeSeconds= $t\n";
my $m = getMemory($cout);
print "Memory=$m\n";

sub readLabel
{
	my ($file, $label) = @_;
	my $value;
	open(FILE, "$file") or die "$0: Cannot open file $file : $!\n";
	while (<FILE>) {
		if (/^$label(.*$)/) {
			$value = $1;
			last;
		}
	}

	close(FILE);

	defined($value) or die "$0: File $file does not contain label $label\n";
	return $value;
}

sub getRunTime
{
	my ($file) = @_;
	my $v1 = readLabel($file, "UnixTimeStart=");
	my $v2 = readLabel($file, "UnixTimeEnd=");
	return ($v2 - $v1);
}

sub getMemory
{
	my ($file) = @_;
	open(FILE, "$file") or die "$0: Cannot open file $file : $!\n";
	my $m = "UNDEFINED";
	while (<FILE>) {
		chomp;
		if (/maximum was (.+$)/) {
			$m = $1;
		}
	}

	close(FILE);
	return $m;
}
