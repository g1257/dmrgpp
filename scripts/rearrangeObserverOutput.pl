#!/usr/bin/perl -w
#
use strict;

my ($file)=@ARGV;

my ($label1,$label2) = ("nupNdown","nUp+nDown"); 

printSuper($file);

rearrange($file,$label1);

rearrange($file,$label2);

sub printSuper
{
	 my ($file)=@_;
	 open(FILE,$file) or die "Cannot open file $file: $!\n";
	 my $found = 0;
	 while(<FILE>) {
	 	if (/superdensity/i) {
			print "###\n";
			print;
			$found  = 1;
			last;
		}
	}
	close(FILE);
	($found) or die "$0: Cannot find superdensity in file $file\n";
}

sub rearrange
{
	my ($file,$label)=@_;
	open(FILE,$file) or die "Cannot open file $file: $!\n";

	while(<FILE>) {
		last if (/^#Using Matrix A:/);
	}
	my $x = $_;
	my $counter = 0;
	while(<FILE>) {
		$x = doOneBlock($label,$x,$counter);
		$counter++;
	}
	close(FILE);

}

sub doOneBlock 
{
	my ($label,$saved,$counter)=@_;
	while(<FILE>) {
		if (/^#/) {
			$saved .= $_;
		}
		if (/^site /) {
			$saved .= $_;
			last;
		}
	}
	my $needsPrinting = 0;
	$needsPrinting = 1 if (/\Q$label/i);

	print $saved if ($needsPrinting and $counter<2); 

	while(<FILE>) {
		last if (/^#/);
		next if (/superdensity/i);
		next if (/^Not found #FERMIONICSIGN in file/);
		next if (/^Ignore prev. error/);
		print if ($needsPrinting);
	}
	return $_;
}

