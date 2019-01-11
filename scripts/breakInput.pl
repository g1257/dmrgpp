#!/usr/bin/perl

use strict;
use warnings;
use strict;

my ($file, $movementsPerFile) = @ARGV;
defined($movementsPerFile) or die "USAGE: $0 inputFile movementsPerFile\n";

open(FILE, "<", "$file") or die "$0: Cannot open $file : $!\n";
my @lines;
while (<FILE>) {
	chomp;
	push @lines, $_;

	# make sure it's a restart run
	if (/SolverOptions/) {
		print STDERR "$0: WARNING: This is not a restart run\n" unless (/restart/);
	}
}

close(FILE);

my @finiteLoops;
# Read all finite loops
my $rest = readFiniteLoops(\@finiteLoops, \@lines);

# Create vector of new finite loops
my @newFl = createNewFiniteLoops(\@finiteLoops, $movementsPerFile);

# output all inputs
my $n = scalar(@newFl);
print STDERR "$0: $n files will be created\n";
my $root = $file;
$root =~ s/\.[^\.]+$//;
for (my $i = 0; $i < $n; ++$i) {
	writeInput($root, $i, $rest, $newFl[$i]);
}

sub writeInput
{
	my ($root, $ind, $rest, $h) = @_;
	my $fout = $root."_$ind".".inp";
	open(FOUT, ">", "$fout") or die "$0: Cannot write to $fout : $!\n";

	print FOUT "$rest\nFiniteLoops 1\n";

	print FOUT $h->{"move"}." ".$h->{"m"}." ".$h->{"option"}."\n";

	close(FOUT);

	print "$0: File $fout has been written\n";
}

sub createNewFiniteLoops
{
	my ($fl, $movementsPerFile) = @_;
	my @newFl;
	my $counter = 0;
	my $hprev;
	my $sum = 0;
	my $n = scalar(@$fl);
	for (my $i = 0; $i < $n; ++$i) {
		my $doBreak = 0;
		my $h = $fl->[$i];
		if ($counter < $movementsPerFile) {
			++$counter;
		} else {
			$doBreak = 1;
		}

		if ($i > 0) {
			if ($h->{"move"} != $hprev->{"move"} ||
			    $h->{"m"} != $hprev->{"m"} ||
			    $h->{"option"} != $hprev->{"option"}) {
				$doBreak = 1;
			}
		}

		if ($doBreak) {
			my $hh = {"move" => $sum, "m" => $hprev->{"m"}, "option" => $hprev->{"option"}};
			push @newFl, $hh;
			$sum = $h->{"move"};
			$counter = 1;
		} else {
			$sum += $h->{"move"};
		}

		$hprev = $h;
	}

	$hprev = $fl->[$n - 1];
	my $hh = {"move" => $sum, "m" => $hprev->{"m"}, "option" => $hprev->{"option"}};
	push @newFl, $hh;
	return @newFl;
}


sub readFiniteLoops
{
	my ($fl, $lines) = @_;
	my ($rest, $blob) = readFiniteLoopsBlob($lines);
	$blob =~ s/^[ \t]+//;
	$blob =~ s/[ \t]+$//;
	my @temp = split/ +/,$blob;
	my $n = scalar(@temp);
	die "$0: FATAL: Too small FiniteLoops line\n" unless $n > 3;
	my $nofl =shift @temp;
	--$n;
	die "$0: FATAL: FiniteLoops: Expected 3n elements, got $n\n" unless ($n % 3 == 0);
	my $nofl2 = int($n/3);
	die "$0: FATAL: Number of Numbers to follow $nofl, expected $nofl2\n" unless ($nofl == $nofl2);

	for (my $i = 0; $i < $nofl; ++$i) {
		my $move = $temp[3*$i];
		my $mvalue = $temp[3*$i + 1];
		my $option = $temp[3*$i + 2];
		my $sign = ($move > 0) ? 1 : -1;
		my $moves = abs($move);
		for (my $j = 0; $j < $moves; ++$j) {
			my $h = {"move" => $sign, "m" => $mvalue, "option" => $option};
			push @$fl, $h;
		}
	}
 
	return ($rest);
}

sub readFiniteLoopsBlob
{
	my ($lines) = @_;
	my $ind = 0;
	my $n = scalar(@lines);
	my $rest = "";
	for (; $ind < $n; ++$ind) {
		last if ($lines[$ind] =~ /^[ \t]*FiniteLoops[ \t]/);
		last if ($lines[$ind] =~ /^[ \t]*FiniteLoops$/);
		$rest .= $lines[$ind]."\n";
	}

	$_ = $lines[$ind++];
	s/^ *FiniteLoops//;
	my $content = $_;
	for (; $ind < $n; ++$ind) {
		last if ($lines[$ind] =~ /[^\d\- \t]/);
		$content .= " ".$lines[$ind];
	}

	for (; $ind < $n; ++$ind) {
		$rest .= $lines[$ind]."\n";
	}

	return ($rest, $content);
}

