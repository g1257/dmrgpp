#!/usr/bin/perl -w

use strict;
use warnings;

my ($observable,$useFlag)=@ARGV;

defined($useFlag) or die "USAGE: $0 observable useFlag\n";

print STDERR "$0: command line: $observable $useFlag\n";

my $nsites=0;
my $flag=0;
my $biggestTimeSeen = 0;
my $obs = "\Q$observable";
my @value;
my ($times,$sites) = (0,0);

while(<STDIN>) {
	chomp;
	next if (/^#/);
	my @temp = split;
	if (scalar(@temp) != 5) { 
		die "$0: Line $. does not have 5 columns\n";
	}

	my ($site,$val,$time,$label,$den) = @temp;

	if (!isAnInteger($site)) {
		die "$0: Line $. site $site is not an integer\n";
	}	
	
	my $timeIndex = int($time*10);
	my $denReal = realPartStrict($den);
	$denReal = 1 if (realPart($denReal) == 0);
	$value[$timeIndex][$site] = realPartStrict($val)/$denReal;
	$times = $timeIndex if ($times < $timeIndex);
	$sites = $site if ($sites < $site);
}

$times++;
$sites++;
print STDERR "$0: Total sites = $sites\n";
print STDERR "$0: Total times = $times\n";

for (my $i = 0; $i < $times; ++$i) {
	my $time = $i*0.1;
	my $c = 0;
	my @temp;
	for (my $j = 0; $j < $sites; ++$j) {
		my $val = $value[$i][$j];
		defined($val) or next;
		$c++;
		$temp[$j] = $val;
	}

	#next if ($c != $sites);

	print "$time ";
	for (my $j = 0; $j < $sites; ++$j) {
		my $val = $temp[$j];
		defined($val) or $val = "0.00";
		print "$val ";
	}

	print "\n";
}
	
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

