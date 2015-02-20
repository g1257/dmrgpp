#!/usr/bin/perl

use strict;
use warnings;
my $counter=0;

my ($site,$label)=@ARGV;
defined($label) or die "USAGE: $0 site label < file\n";

print "#site= $site\n";
print "#label=$label\n";

while (<STDIN>) {
	chomp;
	next unless /\Q$label/;
	if (/^${site} /) {
		my @temp = split;
		my $value = procValue($temp[1]);
		my $time = $temp[2];
		my $superdensity = procValue($temp[4]);
		print "$time  $value $superdensity\n";
	}
}

sub procValue
{
	my ($t)=@_;
	$_=$t;
	s/\(//;
	s/,.*$//;
	return $_;
}

