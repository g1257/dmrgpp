#!/usr/bin/perl

use strict;
use warnings;
use timeObservablesInSitu;


my ($site,$label)=@ARGV;
defined($label) or die "USAGE: $0 site label < file\n";
timeObservablesInSitu::main($site, $label, *STDIN, *STDOUT);

