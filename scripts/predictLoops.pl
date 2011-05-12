#!/usr/bin/perl -w
use strict;
my ($totalT,$deltaT,$sites,$m,$advanceEach);

getInput(\$sites,"Total number of sites","Any",16);
getInput(\$advanceEach,"Advance Each","Any",4);
getInput(\$totalT,"Total Time","Any",2.0);
getInput(\$deltaT,"Delta Time","Any",0.1);
getInput(\$m,"m for finite loops","Any",200);

#print "$totalT $deltaT $sites $m\n";

my $x = $sites/2 - 1;
my $steps = $totalT * $advanceEach;
$steps /=  ($deltaT * $x);
$steps = int($steps) + 1;
my $loops = int($steps/4);

my $buffer="$x $m 0  -$x $m 0  -$x $m 0  $x $m 1\n";
my $count = 4;
for (my $i=0;$i<$loops;$i++) {
	$buffer=$buffer."$x $m 1  -$x $m 1  -$x $m 1  $x $m 1\n";
	$count += 4;
}
print "FiniteLoops $count $buffer\n";

sub getInput
{
	my ($var,$what,$available,$default)=@_;
	print "Please enter the: $what\n";
	print "Available: $available\n";
	print "Default is: $default (press ENTER): ";
	$_=<STDIN>;
	chomp;
	if ($_ eq "" or $_ eq "\n") {
		$_=$default;
	}
	$$var=$_;
}

