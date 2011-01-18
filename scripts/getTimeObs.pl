#!/usr/bin/perl -w

use strict;

my ($site,$file)=@ARGV;

defined($site) or die "$0: Undefined site\n";
defined($file) or die "$0: Undefined file\n";

my $sd = getSuperDensity($site,$file);
open(FILE,$file) or die "$0: Cannot open file $file: $!\n";
while(<FILE>) {
	if (/^site nupNdown\(gs\) nupNdown\(timevector\) time/) {
		last;
	}
}

my $prevT = -1;
while(<FILE>) {
	next if (/^VectorWithOffsets/);
	last if (/^#/);
	my @temp=split;
	if ($temp[0]==$site) {
		my $val = $temp[2];
		$val =~ s/\(//;
		$val =~ s/,.*\)//;
		next if ($prevT == $temp[3]);
		$prevT = $temp[3];
		$val /= $sd;
		print "$temp[3] $val\n";
	}
}
close(FILE);
	
sub getSuperDensity
{
	my ($site,$file)=@_;
	my $sd;
	open(FILE,$file) or die "$0: Cannot open file $file: $!\n";
	while(<FILE>) {
		if (/SuperDensity.*=\(([^,]+),/) {
			$sd = $1;
			last;
		}
	}
	close(FILE);
	defined $sd or die "SuperDensity is not defined\n";
	return $sd;
}
