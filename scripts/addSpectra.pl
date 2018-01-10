#!/usr/bin/perl

use strict;
use warnings;
use utf8;
use Math::Trig;
use OmegaUtils;

use Getopt::Long qw(:config no_ignore_case);

die "$0 is no longer needed; use nFromAkw.pl instead\n";

my $usage = "-f dollarizedInput [-p] [-z]\n";

my $templateInput;
my $isPeriodic;
my $zeroAtCenter = 0;

GetOptions('f=s' => \$templateInput,
           'p' => \$isPeriodic,
           'z' => \$zeroAtCenter) or die "$usage\n";

(defined($templateInput) and defined($isPeriodic)) or die "$0: USAGE: $usage\n";

my $geometry;
my $hptr = {"GeometryKind" => \$geometry};
OmegaUtils::getLabels($hptr,$templateInput);

my $numberKs;
my %a;
foreach my $file (@ARGV) {
	print STDERR "Adding $file\n";
	addToSpectrum(\%a, \$numberKs, $file);
}

#printSpectrum($outSpectrum, \%a);
OmegaUtils::printGnuplot(\%a, $geometry,  $isPeriodic, $zeroAtCenter);

sub printSpectrum
{
	my ($ptr) = @_;
	foreach my $omega (sort {$a <=> $b} keys %$ptr) {
		my $aptr = $ptr->{$omega};
		my @a = @$aptr;
		print "@a \n";
	}
}

sub addToSpectrum
{
	my ($ptr, $ptrN, $file) = @_;

	my $isGood = 1;
	open(FILE, "<", $file) or die "$0: Cannot open $file : $!\n";
	while (<FILE>) {
		next if (/^#/);
		my @temp = split;
		my $n = scalar(@temp);
		if ($n < 2) {
			print STDERR "$0: Ignored line $. in $file, less than 2 cols\n";
			next;
		}

		$$ptrN = $n - 1 if (!defined($$ptrN));
		if ($n - 1 != $$ptrN) {
			$isGood = 0;
			print STDERR "$0: Line $. in $file with wrong cols\n";
			last;
		}

		my $omega = $temp[0];
		my $oldVal = $ptr->{$omega};
		if (defined($oldVal)) {
			for (my $i = 1; $i < $$ptrN; ++$i) {
				$temp[$i] += $oldVal->[$i];
			}
		}

		$ptr->{$omega} = \@temp;
	}

	close(FILE);

	return if ($isGood);

	die "$0: $file with at least 1 line with wrong number of cols, ".$$ptrN." expected\n";
}

