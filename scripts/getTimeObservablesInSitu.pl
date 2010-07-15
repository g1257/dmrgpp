#!/usr/bin/perl -w

use strict;
my $counter=0;

my ($site,$GlobalTimeSteps,$GlobalTau)=@ARGV;

print "#site= $site\n";
print "#timeSteps=$GlobalTimeSteps\n";
print "#tau=$GlobalTau\n";

while (<STDIN>) {
	chomp;
	next unless /P\d\|A\|P\d/;
	if (/^${site} /) {
		$counter = procLine($_,$counter);
	}
}

sub procLine
{
	my ($t,$counter)=@_;
	
	my @temp = split/ +/,$t;
	#print STDERR "Procing temp1=$temp[1]* temp[2]=$temp[2]*\n";
	my $value = procValue($temp[1]);
	my $time = procCounter($temp[2],\$counter);
	#print STDERR "Procing time=$time* value=$value*\n";
	print "$time  $value\n";
	return $counter;
}

sub procCounter
{
	my ($t,$counterPtr)=@_;
	my $x = $$counterPtr;
	$$counterPtr++;
	$$counterPtr = 0 if ($$counterPtr >= $GlobalTimeSteps);
	#print STDERR "procCounter t=$t* GlobalTau=$GlobalTau*\n";
	#print STDERR "c=".$$counterPtr."*\n";
	return $t+($x)*$GlobalTau;	
}

sub procValue
{
	my ($t)=@_;
	$_=$t;
	s/\(//;
	s/,.*$//;
	return $_;
}

