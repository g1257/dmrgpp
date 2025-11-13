#!/usr/bin/perl -w
use warnings;
use strict;

package FunctionParsing;

my $multiLineComment = 0;
my $buffer="";
my $pristineBuffer="";

sub getNextFunction
{
	my ($ffunc,$fh,$doNotPrint)=@_;
	while(<$fh>) {

		my $savedLine = $_;

		s/\/\*.*\*\///;
		s/\/\/.*$//;
		$multiLineComment = 1 if (/^[\t ]*\/\*/);

		if ($multiLineComment) {
			print $savedLine unless ($doNotPrint);
			$multiLineComment = 0 if (/\*\//);
			next;
		}

		# Empty lines:
		my $line = $_;
		$line=~s/[ \t]+//;
		$line=~s/\n//;
		$pristineBuffer .= $savedLine;
		if ($line eq "") {
			print $pristineBuffer unless ($doNotPrint);
			$pristineBuffer="";
			next;
		}

		# Buffer
		$buffer .= $_;

		if (/\{/) {
			my $needsPrinting = procBuffer(\$ffunc,$buffer,$doNotPrint);
			print $pristineBuffer if ($needsPrinting and !$doNotPrint);
			$buffer="";
			$pristineBuffer="";
			next if ($needsPrinting);
			next unless ($doNotPrint);
			my %ggg = %$ffunc;
# 			print STDERR "******".$ggg{"name"}."\n";
			return if ($doNotPrint);
		}

		if (/\}[ \t]*$/ ||  /\:[ \t]*$/ || /\;[ \t]*$/) {
			$buffer = "" ;
			print $pristineBuffer unless ($doNotPrint);
			$pristineBuffer="";
		}
	}
}

sub procBuffer
{
	my ($ffunc,$t,$doNotPrint)=@_;
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
	
	my %func = %$$ffunc;

	$func{"post-qualifier"} = getPostQualifier($t);

	getFirstPart(\%func,$first);
	%$$ffunc = %func;

#   my $name = $func{"name"};
# 	return 1 if (isInVector(\@keywords,$name));
# 
# 	printFunc(\%func);
# 
# 	my @funcArgs=getArgs($second);
# 	for (my $i=0;$i<=$#funcArgs;$i++) {
# 		print STDERR " ~ $funcArgs[$i] ";
# 	}
# 	my $length=computeLength(\%func,\@funcArgs);
# 	print STDERR " LENGTH=$length\n";
# 	print STDERR "\n-----------------------------\n";
# 	return 1 if ($length<$GlobalMaxLine || $#funcArgs<1);
# 	rewriteSig(\%func,\@funcArgs) unless ($doNotPrint);
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

1;
