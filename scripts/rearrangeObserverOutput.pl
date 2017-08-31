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
	 open(FILE, "<", $file) or die "Cannot open file $file: $!\n";
	 my $found = 0;
	 my $saved = "NOT_FOUND";
	 while(<FILE>) {
		if (/superdensity/i) {
			my $newSuperDensity = $_;
			last if (!($saved eq "NOT_FOUND") and
				diffBetweenSuperDensities($newSuperDensity,$saved)<1e-5);
			$saved = $_;
			$found = 1;
		}
	}
	close(FILE);
	($found) or die "$0: Cannot find superdensity in file $file\n";

	print "###\n";
	print $saved;
}

sub diffBetweenSuperDensities
{
	my ($t1,$t2)=@_;
	my $t1N = getNumericSuperDensity($t1);
	my $t2N = getNumericSuperDensity($t2);
	return abs($t2N - $t1N);
}

#SuperDensity(Weight of the timeVector)=(1.23456,0)
sub getNumericSuperDensity
{
	my ($t) = @_;
	my $sLabel = "SuperDensity(Weight of the timeVector)=(";
	$t=~s/\Q$sLabel//;
	$t=~s/,.*$//;
	return $t;
}

sub rearrange
{
	my ($file,$label)=@_;
	open(FILE, "<", $file) or die "Cannot open file $file: $!\n";

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
		if (/superdensity/i) {
			print;
			next;
		}
		next if (/^Not found #FERMIONICSIGN in file/);
		next if (/^Ignore prev. error/);
		print if ($needsPrinting);
	}
	return $_;
}

