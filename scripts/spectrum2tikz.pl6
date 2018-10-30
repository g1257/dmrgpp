#!/usr/bin/perl6

use v6;

sub MAIN($file, $lx)
{
	my $input = open($file, :r);
	my $counter = 0;
	my $tomegas = 10;
	my ($x, $y) = (0, 0);
	while (my $line = $input.get) {
		++$counter;
		my @temp = split(/\s/, $line);
		print "$x $y " ~ @temp[3] ~ "\n";
		++$x;
		if ($x == $lx) {
			$x = 0;
			++$y;
			print "\n";
		}

		#printThisLine(@temp, $lx);
		#last if ($counter == 10*32);
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

