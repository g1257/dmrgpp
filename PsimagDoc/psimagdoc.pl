#!/usr/bin/perl -w

use warnings;
use strict;

our $doxydocs;

use lib 'doxygen/perlmod';

use DoxyDocs;

my $startPtr = $doxydocs;

my $foundHash;
my $verbose = 0;

while(<STDIN>) {
	next if (/\\newcommand\{\\ptex/);

	if (/\\ptexPaste\{([^\}]+)\}/) {
		my $x = $1;
		my @temp=split/,/,$x;
		die "$0: Syntax error for \\ptexPaste line $.\n" if ($#temp<0);
		$foundHash=\$x;
		my $ptr = $startPtr->{"classes"};
		if ($temp[0]=~s/(^!)//) {
			$ptr = $startPtr->{"files"};
		}
		findClassOrFile($ptr,$temp[0],\$foundHash);
		(UNIVERSAL::isa( $foundHash, "HASH" )) or die "$0: Not found file or class $x at line $.\n";
		if ($#temp>=1) {
			my $nameOfKind = $temp[1];
			
			my $kind = "variable";
			if ($nameOfKind=~s/\(\)$//) {
				$kind="function";
			}
			
			findKind($foundHash,$nameOfKind,\$foundHash,$kind);
			(UNIVERSAL::isa( $foundHash, "HASH" ))
				or die "$0: Not found function $x\n";

			printHash($foundHash) if ($verbose);
		}
		$x = $temp[0];
		my $substitution = getDetailed($foundHash);
		$substitution =~ s/!PTEX_THISCLASS/$x/g;
		s/\\ptexPaste\{([^\}]+)\}/$substitution/;
	}
	if (/\\ptexReadFile\{([^\}]+)\}/) {
		print STDERR "Reading $1\n" if ($verbose);
		readFile($1);
		next;
	}
	print;
}

sub readFile
{
	my ($file)=@_;
	open(FILE, "<", $file) or die "Cannot open $file: $!\n";
	while(<FILE>) {
		print;
	}
	close(FILE);
}

sub getDetailed
{
	my ($ptr)=@_;
	my $x = $ptr->{"detailed"}->{"doc"};
	my $buffer = "";
	for my $item (@$x) {
		$buffer .=  procTextual($item);
	}
	return $buffer;
}

sub procTextual
{
	my ($ptr)=@_;
	my $type = $ptr->{"type"};
	return "\n\n" if ($type eq "parbreak");
	return "$type" unless ($type eq "text");
	my $content = $ptr->{"content"};
	defined($content) or $content="";
	return $content;
}

sub printHash
{
	my ($ptr)=@_;
	for my $item (keys %$ptr) {
		my $value = $ptr->{$item};
		print STDERR "key=$item value=$value\n";
	}
}

sub findClassOrFile
{
	my ($ptr,$thisClassOrFile,$foundHash)=@_;
	for my $item (@$ptr) {
		if ($item->{"name"} eq $thisClassOrFile) {
			$$foundHash = $item;
			return;
		}
	}
}

sub findKind
{
	my ($ptr,$thisFunc,$foundHash,$kind)=@_;
	my ($lastKind,$lastName)=("","");
	for my $item (keys %$ptr) {
  		print STDERR "key=$item\n" if ($verbose);

		my $x = $ptr->{$item};
		
		if (UNIVERSAL::isa( $x, "HASH" )) {
			findKind($x,$thisFunc,$foundHash,$kind);
			next;
		} elsif (UNIVERSAL::isa( $x, "ARRAY" )) {
			findKindA($x,$thisFunc,$foundHash,$kind);
			next;
		}
		if ($item eq "kind") {
			$lastKind = $x;
			print STDERR "lastKind=$lastKind\n" if ($verbose);
		}
		if ($item eq "name") {
			$lastName = $x;
			print STDERR "lastName $lastName and thisFunc= $thisFunc\n" if ($verbose);
		}
		
		if ($lastName eq $thisFunc) {
			next if (defined($kind) and !($lastKind eq $kind));
			$$foundHash = $ptr;
 			print STDERR "Found for hash $$foundHash\n" if ($verbose);
			return;
		}
	}
}

sub findKindA
{
	my ($ptr,$thisFunc,$foundHash,$kind)=@_;
	for my $item (@$ptr) {
 		print STDERR "findKindA key=$item\n" if ($verbose);
		my $x = $item;
		if (UNIVERSAL::isa( $x, "HASH" )) {
			findKind($x,$thisFunc,$foundHash,$kind);
		} elsif (UNIVERSAL::isa( $x, "ARRAY" )) {
			findKindA($x,$thisFunc,$foundHash,$kind);
		} else {
 			print STDERR "findKindA value=$x\n" if ($verbose);
		}
	}
}

