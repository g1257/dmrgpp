#!/usr/bin/perl -w

use strict;

# [static const virtual] [type] name(arg1,arg2) [const]

my $multiLineComment = 0;
my $buffer="";
my $pristineBuffer="";

my @keywords=("for","if","while");
my $GlobalSizeOfTab = 2;
my $GlobalMaxLine = 80;

while(<STDIN>) {

	my $savedLine = $_;
	if (/^[\t ]*\/\//) {
		print $savedLine;
		next;
	}
	s/\/\*.*\*\///;
	$multiLineComment = 1 if (/^[\t ]*\/\*/);
	
	if ($multiLineComment) {
		print $savedLine;
		$multiLineComment = 0 if (/\*\//);
		next;
	}

	chomp;

	# Empty lines:
	my $line = $_;
	$line=~s/[ \t]+//;
	$pristineBuffer .= $savedLine;
	if ($line eq "") {
		print $pristineBuffer;
		$pristineBuffer="";
		next;
	}

	# Buffer
	$buffer .= " ".$_;

	if (/\{/) {
		my $needsPrinting = procBuffer($buffer);
		print $pristineBuffer if ($needsPrinting);
		$buffer="";
		$pristineBuffer="";
		next;
	}

	if (/\}[ \t]*$/ ||  /\:[ \t]*$/ || /\;[ \t]*$/) {
		$buffer = "" ;
		print $pristineBuffer;
		$pristineBuffer="";
	}
}

print $pristineBuffer;
$pristineBuffer="";

sub procBuffer
{
	my ($t)=@_;
	return 1 if ($t=~/[ \t]+namespace /);
	return 1 if ($t=~/[ \t]+class /);
	return 1 if ($t=~/\#include /);
	# Find first parens
	my $first="";
	my $second="";
	if ($t=~/(^[^\(]+)\(([^\)]*)\)/) {
		$first = $1;
		$second = $2;
	}
	die "$0: Wrong signature at line $. for $t\n" if ($first eq "");

	
	my %func;
	$func{"post-qualifier"} = getPostQualifier($t);

	getFirstPart(\%func,$first);
	my $name = $func{"name"};
	return 1 if (isInVector(\@keywords,$name));

	printFunc(\%func);

	my @funcArgs;
	getArgs(\@funcArgs,$second);
	for (my $i=0;$i<=$#funcArgs;$i++) {
		print STDERR " ~ $funcArgs[$i] ";
	}
	my $length=computeLength(\%func,\@funcArgs);
	print STDERR " LENGTH=$length\n";
	print STDERR "\n-----------------------------\n";
	return 1 if ($length<$GlobalMaxLine || $#funcArgs<1);
	rewriteSig(\%func,\@funcArgs);
	return 0;
}

sub getFirstPart
{
	my ($f,$first)=@_;
	$_=$first;
	$$f{"level"}=getLevel($_);
	my @temp=split;
	my $n = $#temp+1;
	if ($n==1) { # It's a constructor or destructor
		$$f{"pre-qualifier"}="";
		$$f{"type"}="";
		$$f{"name"}=$temp[0];
		return;
	}
	if ($n==2) { # unqualified function or destructor
		if ($temp[1]=~/^\~/) { # destructor
			$$f{"pre-qualifier"}=$temp[0];
			$$f{"type"}="";
			$$f{"name"}=$temp[1];
			return;
		}
		# unqualified function
		$$f{"pre-qualifier"}="";
		$$f{"type"}=$temp[0];
		$$f{"name"}=$temp[1];
		return;
	}

	# $n>2 --> qualified function
	$$f{"pre-qualifier"} = "";
	for (my $i=0;$i<$n-2;$i++) {
		$$f{"pre-qualifier"} .= $temp[$i];
	}
	$$f{"type"}=$temp[$n-2];
	$$f{"name"}=$temp[$n-1];
}

sub getPostQualifier
{
	my ($t)=@_;
	my $ret = "";
	if ($t=~/\)([^\{\)]*)\{[\t ]*$/) {
		$ret = $1;
		$ret =~ s/[ \t]//g;
	}
	return $ret;
}

sub getLevel
{
	my ($t)=@_;
	return scalar(@{[$t =~ /(\t)/g]});
}

sub getArgs
{
	my ($fa,$t)=@_;
	return if ($t eq "");
	$_ = $t;
	@$fa=split/,/;
	my $n = scalar(@$fa);
	for (my $i=0;$i<$n;$i++) {
		$fa->[$i]=~s/^[\t ]+//;
	}
}

sub isInVector
{
	my ($vec,$what)=@_;
	foreach (@$vec) {
		return 1 if ($what eq $_);
	}
	return 0;
}

sub printFunc
{
	my ($f)=@_;
	foreach (keys %$f) {
		my $val = $$f{"$_"};
		print STDERR "$_=$val ";
	}
}

sub computeLength
{
	my ($f,$fa)=@_;
	my $sum = 0;
	my $sp = 0;
	foreach (keys %$f) {
		my $val = $$f{"$_"};
		next if ($val eq "");
		if ($_ eq "level") {
			$sum += $val*$GlobalSizeOfTab;
			next;
		}
		$sum += length($val);		
		$sp++; # leave a space
	}
	$sp-- if ($sp>0); # substract last space if any
	
	$sum += 2; # add both parens
	
	my $cma = 0;
	foreach (@$fa) {
		next if ($_ eq "");
		$sum += length($_);
		$cma++; # comma between args
	}
	$cma-- if ($cma>0); # substract last comma if any
	return $sum + $cma + $sp;
}

sub rewriteSig
{
	my ($f,$fa)=@_;
	 #[pre-qualifier] [type] name(arg1,
     #                            arg2,
     #                            ...) [post-qualifier]
     #{
	my $level = $$f{"level"};
	printChars($level,"\t");
	$_ = $$f{"pre-qualifier"};
	print "$_ " unless ($_ eq "");
	my $count = length($_);
	$_ = $$f{"type"};
	print "$_ " unless ($_ eq "");
	$count += length($_);
	$_ = $$f{"name"}."(";
	print "$_" unless ($_ eq "");
	$count += length($_);
	$count++;
	# print arguments, not that we have at least 2 arg.
	my $n = scalar(@$fa);
	($n>1) or die "rewriteSig should not have been called, n=$n\n";
	#first arg.
	$_=$fa->[0].",\n";
	print "$_";
	for (my $i=1;$i<$n-1;$i++) {
		printChars($level,"\t");
		printChars($count," ");
		$_=$fa->[$i].",";
		print "$_ ";
	}
	# last arg
	printChars($level,"\t");
	printChars($count," ");
	$_=$fa->[$n-1];
	print "$_".") ".$$f{"post-qualifier"}."\n";
	printChars($level,"\t");
	print "{\n";
	
}

sub printChars
{
	my ($n,$c)=@_;
	for (my $i=0;$i<$n;$i++) { print "$c"; } 
}