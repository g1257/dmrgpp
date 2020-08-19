#!/usr/bin/perl

use strict;
use warnings;
use utf8;

my ($file) = @ARGV;
defined($file) or die "USAGE: $0 filename\n";

my %max = maxFromFile($file);

printMax(\%max);

sub printMax
{
	my ($h) = @_;
	
	foreach my $k (sort {$a <=> $b} keys %$h) {
		my $ptr = $h->{$k};
		my ($omega, $value) = @$ptr;
		print "$k $omega\n";
	}
}

sub maxFromFile
{
	my ($file) = @_;
	my %max;
	open(FILE, "<", "$file") or die "$0: Cannot open $file : $!\n";
	while (<FILE>) {
		my @temp = split;
		next if (scalar(@temp) != 3);
		my $k = $temp[0];
		my $omega = $temp[1];
		my $value = $temp[2];
		
		my $ptr = $max{$k};
		my @newptr = ($omega, $value);
		if (defined($ptr)) {
			my ($oldOmega, $oldValue) = @$ptr;
			$max{$k} = \@newptr if ($value > $oldValue);
		} else {
			$max{$k} = \@newptr;
		}
	}
	
	close(FILE);
	
	return %max;
}
	
