#!/usr/bin/perl

use strict;
use warnings;
use Math::Complex;

my ($beta,$label)=@ARGV;
defined($label) or die "USAGE: $0 beta label < observerOutput.txt\n";

my $flag=0;
my @matrix;
my @sum;
my $counter=0;

while (<STDIN>) {
	if (/^\#Sites=0 1 2/) {
	       $flag=1;
	       next;
	}

	if (/^\#Time=$beta/ and $flag==1) {
		$flag++;
		next;
	}

	if (/^\Q$label/ and $flag==2) {
		matrixRead(\@matrix);
		matrixAdd(\@sum,\@matrix);
		$counter++;
	}

	$flag=0;
}

print STDERR "#Counter=$counter\n";
($counter>0) or die "$0: counter==0\n";

matrixDivide(\@sum,$counter);
print "$label\n";
matrixPrint(\@sum);

sub matrixRead
{
	my ($matrix)=@_;
	$_=<STDIN>;
	my @temp=split;
	my $rows=$temp[0];
	my $cols=$temp[1];
	my $i=0;
	while (<STDIN>) {
		my @temp=split;
		for (my $j=0;$j<scalar(@temp);$j++) {
			my $index = $i+$j*$rows;
			$matrix->[$index]=$temp[$j];
		}

		$i++;
		last if ($i==$rows);
	}
}

sub matrixAdd
{
	my ($matrix1,$matrix2)=@_;
	my $cols = sqrt(scalar(@$matrix2));
	my $rows = $cols;
	for (my $i=0;$i<$rows;$i++) {
		for (my $j=0;$j<$cols;$j++) {
			my $index = $i+$j*$rows;
			if (!defined($matrix1->[$index])) {
				$matrix1->[$index] = Math::Complex->make($matrix2->[$index]);
			} else {
				$matrix1->[$index] += Math::Complex->make($matrix2->[$index]);
			}
		}
	}
}

sub matrixDivide
{
	my ($matrix,$divisor)=@_;
	my $cols = sqrt(scalar(@$matrix));
	my $rows = $cols;
	defined($rows) or die "$0: rows undefined\n";
	for (my $i=0;$i<$rows;$i++) {
		for (my $j=0;$j<$cols;$j++) {
			my $index = $i+$j*$rows;
			defined($matrix->[$index]) or die "$0: matrix $i $j undefined\n";
			$matrix->[$index] /= $divisor;
		}
	}
}

sub matrixPrint
{
	my ($matrix)=@_;
	my $cols = sqrt(scalar(@$matrix));
	my $rows = $cols;
	print "$rows $cols\n";
	for (my $i=0;$i<$rows;$i++) {
		for (my $j=0;$j<$cols;$j++) {
			my $index = $i+$j*$rows;
			print "$matrix->[$index] ";
		}

		print "\n";
	}
}

sub toReal
{
	my ($t) = @_;
	$t=~s/,.*$//;
	$t=~s/\(//;
	return $t;
}

