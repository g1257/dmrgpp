#!/usr/bin/perl -w

use strict;
use warnings;

my ($observable,$useFlag)=@ARGV;

defined($observable) or die "USAGE: $0 observable useFlag\n";
$useFlag=1 if (!defined($useFlag));

my $total = 0;
my @counter;
my @values;
my $nsites=0;
my $flag=0;
my @times;

while(<STDIN>) {
	$flag = 0  if (/ALL LINKS CLEARED/);
	$flag = 1 if (/ALL LINKS SEEN/);
	next unless (/P0/);
	next unless (/$observable/);
	next if ($flag==0 and $useFlag);

	my @temp=split;
	my $denominator = realPart($temp[4]);
	next if ($denominator==0);
	next if (!isAnInteger($temp[0]));
	my $val = realPart($temp[1])/$denominator;
	my $c = $counter[$temp[0]];
	$c=0 if (!defined($c));
	$times[$c]=$temp[2];
	$values[$temp[0]][$c]=$val;
	$counter[$temp[0]]=$c+1;
	$total = $c+1 if ($total<$c+1);
	$nsites = $temp[0] if ($nsites<$temp[0]);
}
$nsites++;

for (my $c=0;$c<$total;$c++) {
	print "$times[$c] ";
	for (my $site=0;$site<$nsites;$site++) {
		my $val = $values[$site][$c];
		$val = 0 if (!defined($val));
		print "$val ";
	}
	print "\n";
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

sub isAnInteger
{
	my ($t)=@_;
	return 0 if (!defined($t));
	return 1 if ($t=~/^[0-9]+$/);
	return 0;
}

