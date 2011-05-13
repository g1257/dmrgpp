#!/usr/bin/perl -w

use strict;
use GetTimeObs;

my ($site,$file,$whatState,$whatObservable)=@ARGV;
my $fh = *STDOUT;
GetTimeObs::main($fh,$site,$file,$whatState,$whatObservable);

