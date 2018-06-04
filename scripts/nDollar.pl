#!/usr/bin/perl

use strict;
use warnings;
use utf8;

use lib ".";
use Ndollar;

my ($file) = @ARGV;
defined($file) or die "USAGE: $0 file\n";

Ndollar::main($file);

