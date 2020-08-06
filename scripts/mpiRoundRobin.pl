#!/usr/bin/perl

use strict;
use warnings;
use utf8;

my ($file, $mpiJobs) = @ARGV;
defined($file) or die "USAGE: $0 filename\n";

my %h;
my $firstNode;
open(FILE, "<", "$file") or die "$0: Cannot open $file : $!\n";
while (<FILE>) {
	chomp;
	next if (/^#/ or $_ eq "");
	my $f = $_;
	$firstNode = $f;
	if (!defined($h{"$f"})) {
		$h{"$f"} = 1;
	} else {
		++$h{"$f"};
	}
}

close(FILE);

my $nodes = scalar(keys %h);
die "$0: Nodes $nodes must be a multiple of mpiJobs $mpiJobs\n" if ($mpiJobs % $nodes != 0);

die "$0: No nodes!\n" if ($nodes == 0 or !defined($firstNode));
my $ppn = $h{"$firstNode"};
my $repeat = $mpiJobs/$nodes;

my $firstcall = 1;
for (my $i = 0; $i < $repeat; ++$i) {
	foreach my $key (sort keys %h) {
		my $n = $h{"$key"};
		next if ($n == 0);
		if ($firstcall) {
			$firstcall = 0;
		} else {
			print ",";
		}
	
		print "$key";
	}
}


