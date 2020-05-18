#!/usr/bin/perl
# All credit goes to Nirav P.
# All errors go to G.A.

use strict;
use warnings;
use utf8;
use lib ".";
use Honeycomb;

my ($templateInput, $observable) = @ARGV;
defined($observable) or die "USAGE: $0 templateInput observable\n";

my $honey= Honeycomb::init($templateInput);

defined($honey) or die "$0: No honey\n";

my $n1neigh = $honey->{"n1neigh"};
my $n = $honey->{"info"}->{"n"};

if ($observable eq "je") {
	my $str = getJe($n1neigh, $n);
	print "je=$str\n";
	exit(0);
} elsif ($observable eq "js") {
	my $str = getJs($n1neigh, $n);
	print "je=$str\n";
	exit(0);
} else {
	die "$0: Don't know how to print observable $observable\n";
}

sub getJe
{
	my ($n1neigh, $n) = @_;
	# Let $i be a site in the honeycomb lattice, then
	# what is the neighbor of $i in direction $dir, it's
	# $j = n1neigh->[$i + $dir*$n]
	# where $dir = 0, 1, or 2.
	my $sites = scalar(@$n1neigh);
	
	#Just an example for now
	# This is \sum_i sz[i] * sz[i+x] * |gs>
	my $str = "";
	for (my $i = 0; $i < $sites; ++$i) {
		my $j = $n1neigh->[$i + 0*$n]; # neighbor of $i in the x direction $dir = 0
		$str .= "+" if ($i > 0); # no leading plus
		$str .= "sz[$i]*sz[$j]*|gs>";
	}

	return $str;
}

sub getJs
{
	my ($n1neigh, $n) = @_;
	return "spin current NOT DONE YET (sorry) because I'm a slacker :-(\n";
}

