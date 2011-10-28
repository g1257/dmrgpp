#!/usr/bin/perl -w

use strict;

# [static const virtual] [type] name(arg1,arg2) [const]

my $multiLineComment = 0;
my $buffer="";
my $pristineBuffer="";

my @keywords=("for","if","while","struct","class","namespace");
my $GlobalSizeOfTab = 2;
my $GlobalMaxLine = 80;
my @GlobalCvQualifiers=("const","virtual","volatile","static");

while(<STDIN>) {

	my $savedLine = $_;
	
	s/\/\*.*\*\///;
	s/\/\/.*$//;
	$multiLineComment = 1 if (/^[\t ]*\/\*/);
	
	if ($multiLineComment) {
		print $savedLine;
		$multiLineComment = 0 if (/\*\//);
		next;
	}

	# Empty lines:
	my $line = $_;
	$line=~s/[ \t]+//;
	$line=~s/\n//;
	$pristineBuffer .= $savedLine;
	if ($line eq "") {
		print $pristineBuffer;
		$pristineBuffer="";
		next;
	}

	# Buffer
	$buffer .= $_;

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
	return 1 if ($t=~/\#define /);
	return 1 if ($t=~/\#include /);
	# Find first parens
	my $first="";
	my $second="";
	if ($t=~/(^[^\(]+)\(([^\)]*)\)/) {
		$first = $1;
		$second = $2;
	}
	return 1 if ($first eq "");

	
	my %func;
	$func{"post-qualifier"} = getPostQualifier($t);

	getFirstPart(\%func,$first);
	my $name = $func{"name"};
	return 1 if (isInVector(\@keywords,$name));

	printFunc(\%func);

	my @funcArgs=getArgs($second);
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
	$$f{"level"}=getLevel($first);
	$first=~s/\n/\@ /g;
	$_=$first;
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
	$_=$first;
	s/[^ \t]+[\t ]*$//; # kill name
	s/[^ \t]+[\t ]*$//; #kill type
	s/^[\t ]*//;
	s/[\t ]*$//;
	$$f{"pre-qualifier"} = $_;
	$$f{"name"}=$temp[$n-1];
	$$f{"type"}=$temp[$n-2];
}

sub getPostQualifier
{
	my ($t)=@_;
	my $ret = "";
	$t=~s/\n/\@/g;
	if ($t=~/\)(.*$)/) {
		$ret = $1;
		$ret =~ s/[ \t]+$//g;
	}
	return $ret;
}

sub getLevel
{
	my ($t)=@_;
	$_=$t;
	if ($t=~/\n/) {
		my @temp=split/\n/;
		$_=$temp[$#temp];
	}
	#print STDERR "HERE $_\n";
	my $counter=0;
	while(s/^\t//) {
		#print STDERR "*".$1."*\n";
		$counter++;
	}
	
	return $counter;
}

sub getArgs
{
	my ($t)=@_;
	return if ($t eq "");
	$_ = $t;
	my @fa=split/,/;
	my $n = $#fa+1;
	for (my $i=0;$i<$n;$i++) {
		$fa[$i]=~s/^[\t ]+//;
		$fa[$i]=~s/[\t ]+$//;
		$fa[$i]=~s/[\n\r]//g;
	}
	return glueArgs(\@fa);
}

sub glueArgs
{
	my ($fa)=@_;
	my @ret;
	my $buffer="";
	my $openFlag=0;
	my $n = scalar(@$fa);
	my $counter=0;
	for (my $i=0;$i<$n;$i++) {
		$buffer .= "," unless ($buffer eq "");
		$_ = $fa->[$i];
		s/^[\t ]*//;
		s/[\t ]*$//;
		$buffer .= $_;
		$openFlag += countChars($fa->[$i],"<");
		$openFlag -= countChars($fa->[$i],">");
		if ($openFlag==0) {
			$ret[$counter++]=$buffer;
			$buffer="";
		}
	}
# 	$ret[$counter++]=$buffer;
# 	$buffer="";
	return @ret;
}

sub countChars
{
	my ($t,$c)=@_;
	$_=$t;
	my $counter=0;
	while(s/\Q$c//) {
		$counter++;
	}
	return $counter;
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
	$$f{"pre-qualifier"} = templatedFunction($$f{"pre-qualifier"},$level);
	printChars($level,"\t");
	$_ = $$f{"pre-qualifier"};
	print "$_ " unless ($_ eq "");
	my $count = length($_);
	$_ = $$f{"type"};
	print "$_ " unless ($_ eq "");
	$count += length($_);
	$count-- if ($_ eq "");
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
		$_=$fa->[$i].",\n";
		s/^[ \t]+//;
		print "$_";
	}

	# last arg
	printChars($level,"\t");
	printChars($count," ");
	$_=$fa->[$n-1];
	s/^[ \t]+//;
	s/[ \t]*$//;
	print "$_".")";
	printPostQualifier($$f{"post-qualifier"},$level);
	#s/^[ \t]+//;
	#s/[\r\n]//g;
	#print " $_" unless ($_ eq "");
	#print "\n";
	#printChars($level,"\t");
	#print "{\n";
	
}

sub printPostQualifier
{
	my ($pq,$level)=@_;

	return if ($pq eq "");

	if (!($pq=~/\@/)) {
		$pq=~s/[\t ]*$//;
		print "$pq\n";
		return;
	}

	$_=$pq;
	my @temp = split/\@/;
	$_=shift @temp;
	print  "$_\n";
	foreach (@temp) {
		#printChars($level,"\t");
		s/[\t ]*$//;
		print "$_\n"
	}
}

sub templatedFunction
{
	my ($t,$level)=@_;
	
	return $t unless ($t=~/[\t ]*template[ \<]/);
	my ($templ,$cvq)=getCvQualifiers($t);
	printChars($level,"\t");
	print "$templ\n";
	return $cvq;
}

sub getCvQualifiers
{
	my ($t)=@_;
	my ($templ,$cvq)=($t,"");
	$_=$t;
	while(s/([^ \t]+)[\t ]*$//) {
		if (isInVector(\@GlobalCvQualifiers,$1)) {
			$cvq=$1." ";
			$templ=$_;
		} else {
			last;
		}
	}
	$cvq=~s/[\t ]*$//;
	$templ=~s/\@/\n/g;
	$templ=~s/[\n \t]*$//;
	return ($templ,$cvq);
}

sub printChars
{
	my ($n,$c)=@_;
	for (my $i=0;$i<$n;$i++) { print "$c"; } 
}
