#!/usr/bin/perl

use strict;
use warnings;
use utf8;

my @extensions=("*.hh", "*.cc", ".cpp");

my $opts = " --style=file -i ";
foreach my $ext (@extensions) {
	my $cmd = "find . -iname \"$ext\" -exec clang-format  $opts '{}' \\;";
	print STDERR "$0: Formated all files with extension $ext\n";
	system($cmd);
}

