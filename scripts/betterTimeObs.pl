#!/usr/bin/perl

use strict;
use warnings;
use utf8;

my ($file, $label) = @ARGV;
defined($label) or die "USAGE: $0 file label\n";

print "#CmdLine: $file $label\n";

my @a;
my %times;
open(FIN, "<", $file) or die "$0: Cannot open $file : $!\n";
while (<FIN>) {
	next unless /\Q$label/;
	next if (/CmdLine/);

	my @tmp = split;
	next if (scalar(@tmp) < 5);

	my $site = $tmp[0];
	my $value = realPartOf($tmp[1]);
	my $time = $tmp[2];
	my $superdensity = realPartOf($tmp[4]);

	$a[$site]{"$time"} = {"value" => $value, "superdensity" => $superdensity};
	$times{"$time"} = 1;
}

close(FIN);

my $nsites = scalar(@a);
print STDERR "$0: Found $nsites sites in $file\n";

for my $time (sort keys %times) {
	for (my $site = 0; $site < $nsites; ++$site) {
		my $pair = $a[$site]{"$time"};
		if (!defined($pair)) {
			print "$time $site -100 -100\n";
			next;
		}

		print "$time $site ".$pair->{"value"}." ".$pair->{"superdensity"}."\n";
	}
}

sub realPartOf
{
	my ($x) = @_;
	$_ = $x;
	s/\,.*$//;
	s/\(//;
	return $_;
}

