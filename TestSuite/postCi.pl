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
	procMemcheck($n);
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

sub procMemcheck
{
	my ($n) = @_;
	my $file1 = "tests/memcheck$n.txt";
	my $mode = "OK";
	my $extra = "UNDEFINED";
	my %lost;
	open(FILE,$file1) or die "$0: Cannot open $file1 : $!\n";
	while (<FILE>) {
		if (/FATAL: You are requesting/ || /mandatory label/) {
			$mode = "FATAL";
			$extra = $_;
			chomp($extra);
			last;
		}

		next unless (/^==/);
		if (/invalid/i) {
			$mode = "invalid";
			last;
		}

		if (/uninitial/i) {
			$mode = "uninitialized";
			last;
		}
		
		if (/([^ ]+) lost: ([^ ]+) bytes/) {
			my $val = $2;
			my $name = $1;
			$val =~ s/,//;
			$lost{"$name"}= $val;
			next;
		}
	}

	close(FILE);

	if ($mode eq "FATAL") {
		print "$0: ATTENTION TEST $n couldn't run because of $extra\n";
		return;
	}
	
	foreach my $key (keys %lost) {
		my $val = $lost{"$key"};
		next if ($val == 0);
		print "$0: ATTENTION TEST $n $key lost ".$lost{"$key"}." bytes\n";
	}

	return if ($mode eq "OK");
	print "$0: ATTENTION: TEST $n has memcheck mode $mode\n";
}

