#!/usr/bin/perl -w
use strict;
use AkimaSumBatch;

my ($root,$total,$ext)=@ARGV;
$ext = ".txt" if !defined($ext);
my ($start,$end,$points)=(0,3.0,100);
my $fh = *STDOUT; # write to stdout
AkimaSumBatch::main($fh,$start,$end,$points,$root,$total,$ext);


