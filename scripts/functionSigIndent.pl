#!/usr/bin/perl -w

use strict;

# [static const] [type] name(arg1,arg2)

my $multiLineComment = 0;
my $buffer="";

my @keywords=("for","if","while");

while(<STDIN>) {

	next if (/^[\t ]*\/\//);
	s/\/\*.*\*\///;
	$multiLineComment = 1 if (/^[\t ]*\/\*/);
	$multiLineComment = 0 if (/\*\//);
	next if ($multiLineComment);
	chomp;

	# Empty lines:
	my $line = $_;
	$line=~s/[ \t]+//;
	next if ($line eq "");
	
	# Buffer
	$buffer .= " ".$_;
	$buffer = "" if (/\;[ \t]*$/);
	if (/\{/) {
		procBuffer($buffer);
		$buffer="";
		next;
	}
	$buffer = "" if (/\}[ \t]*$/);
	$buffer = "" if (/\:[ \t]*$/);
}

sub procBuffer
{
	my ($t)=@_;
	return if ($t=~/[ \t]+namespace /);
	return if ($t=~/[ \t]+class /);
	return if ($t=~/\#include /);
	# Find first parens
	my $first="";
	my $second="";
	if ($t=~/(^[^\(]+)\(([^\)]*)\)/) {
		$first = $1;
		$second = $2;
	}
	die "$0: Wrong signature at line $. for $t\n" if ($first eq "");

	my %func;
	getFirstPart(\%func,$first);
	my $name = $func{"name"};
	return if (isInVector(\@keywords,$name));
	#print STDERR "$first <-------> $second\n";
	print "c-v=".$func{"c-v"}." type=".$func{"type"}." name=".$func{"name"}." ";
	my @funcArgs;
	getArgs(\@funcArgs,$second);
	for (my $i=0;$i<=$#funcArgs;$i++) {
		print " ~ $funcArgs[$i] ";
	}
	print "\n-----------------------------\n";
}

sub getFirstPart
{
	my ($f,$first)=@_;
	$_=$first;
	my @temp=split;
	my $n = $#temp+1;
	if ($n==1) { # It's a constructor or destructor
		$$f{"c-v"}="";
		$$f{"type"}="";
		$$f{"name"}=$temp[0];
		return;
	}
	if ($n==2) { # unqualified function or destructor
		if ($temp[1]=~/^\~/) { # destructor
			$$f{"c-v"}=$temp[0];
			$$f{"type"}="";
			$$f{"name"}=$temp[1];
			return;
		}
		# unqualified function
		$$f{"c-v"}="";
		$$f{"type"}=$temp[0];
		$$f{"name"}=$temp[1];
		return;
	}
	# $n>2 --> qualified function
	$$f{"c-v"} = "";
	for (my $i=0;$i<$n-2;$i++) {
		$$f{"c-v"} .= $temp[$i];
	}
	$$f{"type"}=$temp[$n-2];
	$$f{"name"}=$temp[$n-1];
}

sub getArgs
{
	my ($fa,$t)=@_;
	return if ($t eq "");
	$_ = $t;
	@$fa=split/,/;
}

sub isInVector
{
	my ($vec,$what)=@_;
	foreach (@$vec) {
		return 1 if ($what eq $_);
	}
	return 0;
}

		

	