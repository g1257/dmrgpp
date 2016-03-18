#!/usr/bin/perl

use strict;
use warnings;
use Ci;

my ($min,$max) = @ARGV;

my @tests;
Ci::getTests(\@tests);

my $total = scalar(@tests);
for (my $i = 0; $i < $total; ++$i) {
	my $n = $tests[$i];
	next if (defined($min) and $n < $min);
	next if (defined($max) and $n > $max);
	procTest($n);
}
sub procTest
{
	my ($n) = @_;
	procData($n);
}

sub procData
{
	my ($n) = @_;
	my $file1 = "tests/data$n.txt";
	my $file2 = "oldTests/data$n.txt";
	(-r "$file1") or return;
	(-r "$file2") or return;
	my $cmd = "diff $file1 $file2";
	my @version = ("???","???");
	open(PIPE,"$cmd |") or return;
	while (<PIPE>) {
		if (/([\<\>]) DMRG\+\+ version (.*$)/) {
			my $tmp = ($1 eq "<") ? 0 : 1;
			$version[$tmp] = $2;
			next;
		}
	}

	close(PIPE);

	print "New Version $version[0], Old Version $version[1]\n";
}

