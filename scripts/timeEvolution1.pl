#!/usr/bin/perl -w

use strict;
use warnings;

my ($obs,$useFlag)=@ARGV;

defined($useFlag) or die "USAGE: $0 observable useFlag\n";

print STDERR "$0: command line: $obs $useFlag\n";

my $nsites=0;
my $flag=0;
my $biggestTimeSeen = 0;

while(<STDIN>) {
	my $line = $_;
	$flag = 1 if (/^ALL/);
	next unless (/\Q$obs/);
	next if ($flag==0 and $useFlag==1);

	my @temp=split;
	my $denominator = realPartStrict($temp[4]);
	next if ($denominator==0);
	next if (!isAnInteger($temp[0]));

	print $line;

	my $time = $temp[2];
	my $site = $temp[0];
	if ($biggestTimeSeen > $time) {
		die "$0: Error at $. for time=$time\n$_\n";
	}

	$biggestTimeSeen = $time if ($biggestTimeSeen < $time);

	$nsites = $site if ($nsites < $site);
}

$nsites++;

print STDERR "$0: Largest site seen = $nsites\n";
print STDERR "$0: Largest time seen = $biggestTimeSeen\n";

sub realPartStrict
{
	my ($t)=@_;
	my $it = imagPart($t);
	if (abs($it) > 1e-6) {
		die "$0: $t is has non-zero imaginary part\n";
	}

	return realPart($t);
}

sub realPart
{
	my ($t)=@_;
	$_=$t;
	return "-1" if (!defined($_));
	s/\(//;
	s/\,.*$//;
	return $_;
}

sub imagPart
{
        my ($t)=@_;
        $_=$t;
        return "-1" if (!defined($_));
        s/\)//;
        s/^.*\,//;
        return $_;
}

sub isAnInteger
{
	my ($t)=@_;
	return 0 if (!defined($t));
	return 1 if ($t=~/^[0-9]+$/);
	return 0;
}

