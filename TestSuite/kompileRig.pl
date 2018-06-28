#!/usr/bin/perl

use strict;
use warnings;
use utf8;

use lib ".";
use Combinatorial;

my @levels;
my $counter = 0;
$levels[$counter++] = ["", "CPPFLAGS += -DUSE_FLOAT"];
$levels[$counter++] = ["", "CXX = clang++ -mtune=native"];

Combinatorial::main(\@levels, \&myCallback);

sub myCallback
{
	my ($items) = @_;
	my $n = scalar(@$items);
	for (my $i = 0; $i < $n; ++$i) {
		print "$items->[$i]\n";
	}

	print "---------------------\n";
}
