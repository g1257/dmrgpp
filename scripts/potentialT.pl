#!/usr/bin/perl -w
use strict;
use warnings;

print "Enter number of sites\n";
print "Available: Any\n";
print "Default is 16 (press ENTER) ";
$_=<STDIN>;
chomp;
if ($_ eq "") {
	$_= 16;
}       
my $sites = $_;

print "For E_0*i*cos(omega*time)\n";
print "Enter E_0\n";
print "Available: Any\n";
print "Default is -0.2 (press ENTER) ";
$_=<STDIN>;
chomp;
if ($_ eq "") {
	$_= -0.2;
}
my $E0 = $_;

print "Enter omega\n";
print "Available: Any\n";
print "Default is 0.8 (press ENTER) ";
$_=<STDIN>;
chomp;
if ($_ eq "") {
	$_= 0.8;
}       
my $omega = $_;

print "Copy this into your input file in the correct place:\n\n";

print "potentialT $sites ";
for (my $i=0;$i<$sites;$i++) {
	my $val = $E0*($i+1);
	print "$val ";
}
print "\n";
print "omega=$omega\n";




