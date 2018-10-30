#!/usr/bin/perl6

use v6;

sub MAIN($file, $lx)
{
	my $input = open($file, :r);
	my ($x, $y) = (0, 0);
	while (my $line = $input.get) {
		my @temp = split(/\s/, $line);
		#print "$x $y " ~ @temp[3] ~ "\n";
		printThisLine(@temp, $lx);
		++$x;
		if ($x == $lx) {
			$x = 0;
			++$y;
			print "\n";
		}

	}

	$input.close;
}

sub printThisLine($a, $lx)
{
	my $kx = $a.[0];
	my $omega = $a.[1];
	my $imag = $a.[3];
	print "$kx $omega $imag\n";
}

