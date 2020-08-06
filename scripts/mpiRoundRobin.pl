#!/usr/bin/perl

use strict;
use warnings;
use utf8;

my ($file) = @ARGV;
defined($file) or die "USAGE: $0 filename\n";

my %h;
open(FILE, "<", "$file") or die "$0: Cannot open $file : $!\n";
while (<FILE>) {
	chomp;
	next if (/^#/ or $_ eq "");
	my $f = $_;
	if (!defined($h{"$f"})) {
		$h{"$f"} = 1;
	} else {
		++$h{"$f"};
	}
}

close(FILE);

my $firstcall = 1;
my $maxn = 1;
while ($maxn > 0) {
	$maxn = 0;
	foreach my $key (sort keys %h) {
		my $n = $h{"$key"};
		next if ($n == 0);
		$maxn += $n;
		if ($firstcall) {
			$firstcall = 0;
		} else {
			print ",";
		}
	
		print "$key";
		--$h{"$key"};
	}
}


