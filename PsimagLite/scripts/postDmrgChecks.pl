#!/usr/bin/perl

use strict;
use warnings;
use utf8;

# BEGIN: values to check go here
my $test2 = {n => 2, 'Found lowest eigenvalue' => -9.51754};
my $test3 = {n => 3, 'Found lowest eigenvalue' => -21.5102};
# END (values to check)

my ($dir) = @_;
defined($dir) or $dir = ".";
my ($EC_UNDEFINED, $EC_NUMERIC, $EC_STRING) = (-1, -2, -3);

my @tests = ($test2, $test3);

my $flag = 0;
for my $test (@tests) {
	my $log = checkLabelAndValue($test, $dir);
	if ($log ne "") {
		print STDERR "$log\n";
		++$flag;
	}
}

exit($flag);

sub checkLabelAndValue
{
	my ($gold, $dir) = @_;
	my $n = $gold->{n}; # test number
	my $eps = $gold->{eps};
	$eps = 1e-5 if !defined($eps);
	my $file = "runForinput$n.cout"; # construct filename
	my @labels = keys %$gold; # all labels
	my %steel = findValues($dir, $file, \@labels); # steel values
	my $str = "";
	foreach my $key (@labels) {
		next if ($key eq "n");
		next if ($key eq "eps");
		my $val = $gold->{$key}; # gold value
		my $diff = findDifference($steel{$key}, $val);
		my $can_compare = $diff->{can_compare};
		my $difference = $diff->{difference};
		$diff->{label} = $key;
		$diff->{n} = $n;
		if (!$can_compare || $difference > $eps) {
			$str .= logFailure($diff, $steel{$key}, $val);
		}
	}

	return $str;
}

sub logFailure
{
	my ($diff, $steel, $gold) = @_;
	my $n = $diff->{n};
	my $label = $diff->{label};
	my $str = "$0: [$n] failed for label $label: ";
	my $tmp = findExplanation($diff, $steel, $gold);
	return "$str $tmp\n";
}

sub findExplanation
{
	my ($diff, $steel, $gold) = @_;
	my $tmp = "";
	if (!$diff->{can_compare}) {
		$tmp = "Cannot compare.";
	} else {
		$tmp = "Difference ".$diff->{difference}." is too big.";
	}

	my $gld = ($gold) ? $gold : "undefined";
	my $stl = ($steel) ? $steel : "undefined";
	$tmp .= " gold = $gld, steel = $stl\n";

	return $tmp;
}

sub findDifference
{
	my ($a, $b) = @_;
	my $isNa = isNumeric($a);
	my $isNb = isNumeric($b);
	my $can_compare = 0;
	my $diff = 0;
	if ($isNa == $isNb and $isNa == $EC_NUMERIC) {
		$diff = abs($a - $b);
		$can_compare = 1;
	}

	return {isA => $isNa, isB => $isNb, can_compare => $can_compare, difference => $diff};
}

sub isNumeric
{
	my $a = shift;
	defined($a) or return $EC_UNDEFINED;
	return ($a =~ /^[\d\.eE\-\+]+$/) ? $EC_NUMERIC : $EC_STRING;
}

# I call "steel" the values produced by this run
# as opposed to the reference or correct or gold values
# stored by this script
sub findValues
{
	my ($dir, $filename, $labels) = @_;
	my %steel;
	my $file = "$dir/$filename";
	open(my $fh, "<", $file) or die "$0: Cannot open $file : $!\n";
	while (<$fh>) {
		foreach my $label (@$labels) {
			# It will overwrite so that
			# only last value found is stored
			if (/$label *\= *([^ ]+)/) {
				$steel{$label} = $1;
			}
		}
	}

	close($fh);
	return %steel;
}

