#!/usr/bin/perl

use strict;
use warnings;
use utf8;
use lib ".";
use MinCompile;

my @models = @ARGV;

if (scalar(@models) == 0) {
	my $usage = "USAGE: $0 what [word0 word1 ...]\n";
	$usage .= "\t what can be models or targets\n";
	$usage .= "\t without words it restores file\n";
	die($usage);
}

my $what = shift @models;

MinCompile::main($what, \@models);

