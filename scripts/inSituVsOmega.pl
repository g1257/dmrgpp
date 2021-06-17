#!/usr/bin/perl

=pod
USAGE is inSituVsOmega.pl runForinput label [total]
OUTPUT is 2 columns
Omega Values
=cut

use warnings;
use strict;
use utf8;
use MIME::Base64;

my ($root, $label, $total, $wantsMax) = @ARGV;
defined($label) or die "USAGE: $0 runForinput label [total] [wantsMax]\n";

#Avoid infinite loops here:
my $hasTotal = 1;
if (!defined($total)) {
	$hasTotal = 0;
	$total = 1000;
}

#\d+.cout
my $start = 0;
if ($root =~ /(\d+)\.cout/) {
	$start = $1;
	$root =~ s/\d+\.cout//;
}

print STDERR "#$root $label ";
print STDERR "$total" if ($hasTotal);
print STDERR "\n";
print STDERR "#Note value may take two columns if it is complex\n";
print STDERR "#Omega Value\n";
for (my $i = $start; $i < $total; ++$i) {
	my $file = $root."$i".".cout";
	if (-r "$file") {
		my $input = loadInputFromCout($file);
		my $omega = extractOmegaFromInput($input);
		my ($value, $max) = extractValueFromCout($file, $label);
		$max = "" unless ($wantsMax);
		print "$omega $value $max\n";
		next;
	}

	# File not readable or doesn't exist
	last unless ($hasTotal);
}

sub extractValueFromCout
{
	my ($file, $label) = @_;
	open(FILE, "<", $file) or die "$0: Cannot open $file : $!\n";
	my $value;
	my $prev;
	my $max = 0;
	while (<FILE>) {
		chomp;
		if (/$label/) {
			my @temp = split;
			my $n = scalar(@temp);
			next if ($n != 5);
			$value = $temp[1];
			$value =~s/\(//;
			$value =~s/\)//;
			$value =~s/,/ /;
		}

		if (defined($prev)) {
			my $diff = findAbsDiff($prev, $value);
			$max = $diff if ($diff > $max);
		}

		$prev = $value;
	}

	close(FILE);

	defined($value) or die "$0: $label not found in $file\n";
	return ($value, $max);
}

sub loadInputFromCout
{
	my ($file) = @_;
	open(FILE, "<", $file) or die "$0: Cannot open $file : $!\n";

	my $input;
	while (<FILE>) {

		if (/PsiApp::echoBase64: Echo of ([^ ]+) /) {
			$_ = <FILE>;
			chomp;
			$input = decode_base64($_);
			last;
		}
	}

	close(FILE);

	defined($input) or die "$0: PsiApp::echoBase64: not found in $file\n";
	return $input;
}

sub extractOmegaFromInput
{
	my ($input) = @_;
	my @lines = split/\n/, $input;
	foreach my $line (@lines) {
		if ($line=~/CorrectionVectorOmega=(.+)/) {
			return $1;
		}
	}

	die "$0: Did not find CorrectionVectorOmega=\n";
}

sub findAbsDiff
{
	my ($a, $b) = @_;
	my @temp1 = split/ /, $a;
	my @temp2 = split/ /, $b;
	my $n = scalar(@temp1);
	(scalar(@temp2) == $n) or die "$0: findAbsDiff different sizes for $a and $b\n";
	($n > 0 and $n < 3) or die "$0: findAbsDiff: sizes not 1 or 2 for $a and $b\n";
	return abs($a - $b) if ($n == 1);
	my $x = $temp1[0] - $temp2[0];
	my $y = $temp1[1] - $temp2[1];
	return sqrt($x*$x + $y*$y);
}

