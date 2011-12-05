#!/usr/bin/perl -w

use warnings;
use strict;

=pod
Usage:
perl ptex.pl < source.ptex > destination.tex
=cut

my $ptexStartKeyword = "!PTEX-START";
my $ptexEndKeyword = "!PTEX-END";
my $ptexThisClassKeyword="!PTEX_THISCLASS";

my %GlobalMacros=();
my $GlobalPtexOpen=0;
my $GlobalLabel="";
my $GlobalBuffer="";

loadMacros();

replaceMacros();

sub replaceMacros
{
	while(<STDIN>) {
		if (/\\pasteRaw\{([^\}]+)\}/) {
			die "Undefined macro for label $1\n" if (!defined($GlobalMacros{$1}));
		}
		s/\\pasteRaw\{([^\}]+)\}/$GlobalMacros{$1}/;
		next if (/\!PTEX\-/);
		print;
	}
}

sub loadMacros
{
	open(PIPE,"find ../src -iname \"*.h\" |") or die "Cannot open pipe: $!\n";
	while(<PIPE>) {
		chomp;
		procThisFile($_);
	}
	close(PIPE);
}

sub procThisFile
{
	my ($file)=@_;
	$GlobalPtexOpen=0;
	$GlobalLabel="";
	$GlobalBuffer="";
	open(FILE,$file) or die "Cannot open file $file: $!\n";
	while(<FILE>) {
		procThisLine($_,$file);
	}
	close(FILE);
}

sub procThisLine
{
	my ($line,$file)=@_;
	my $class = $file;
	$class=~s/\.h$//;
	if ($GlobalPtexOpen) {
		$line=~s/$ptexThisClassKeyword/$class/g;
		if ($line=~/$ptexEndKeyword/) {
			$GlobalPtexOpen=0;
			$GlobalMacros{"$GlobalLabel"}=$GlobalBuffer;
# 			print STDERR "Macro found for $GlobalLabel\n";
			$GlobalBuffer="";
			return;
		}
		$GlobalBuffer .= $line;
		return;
	}

	if ($line=~/$ptexStartKeyword +(.*$)/) {
		$GlobalLabel = $1;
# 		print STDERR "GL=$1\n";
		$GlobalLabel=~s/ //g;
		die "Empty level on line $.\n" if ($GlobalLabel eq "");
		$GlobalPtexOpen=1;
		$GlobalBuffer="";
	}
}


