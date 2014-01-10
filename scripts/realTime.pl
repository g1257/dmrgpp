#!/usr/bin/perl -w
use strict;

my ($label)=@ARGV;
my ($initial,$final);

$label = "\ 2014\$" if (!defined($label));

while(<STDIN>) {
	if (/$label/ and !defined($initial)) {
		$initial = $_;
		next;
	}
	if (/$label/) {
		$final = $_;
		last;
	}
}

die "Undefined initial\n" if (!defined($initial));
die "Undefined final\n" if (!defined($final));


my $rt = getRealTime($initial,$final);
print "$rt\n";

sub getRealTime
{
	my ($i,$j)=@_;
	my $ti = getTimeInSeconds($i);
	my $tj = getTimeInSeconds($j);
	return $tj-$ti;
}

sub getTimeInSeconds
{
	my ($t)=@_;
	#Fri Nov 13 20:36:05 2009
	my @temp = split(/ +/,$t);
	die "Error in $t\n" if ($#temp<3);
	return timeInSeconds($temp[3])+$temp[2]*86400;
}

sub timeInSeconds
{
	#20:36:05
	my ($t)=@_;
	my @temp = split(/:/,$t);
	die "Error in $t (timeInSeconds)\n" if ($#temp<2);
	return $temp[0]*3600 + $temp[1]*60 + $temp[2];
}


