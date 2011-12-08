#!/usr/bin/perl -w

use warnings;
use strict;

our $doxydocs;

use DoxyDocs;

my $startPtr = $doxydocs;

my $foundHash;

while(<STDIN>) {
	next if (/\\newcommand\{\\ptex/);

	if (/\\ptexPaste\{([^\}]+)\}/) {
		my $x = $1;
		my @temp=split/,/,$x;
		die "$0: Syntax error for \\ptexPaste line $.\n" if ($#temp<0);
		my $ptr = $startPtr->{"classes"};
		$foundHash=\$x;
		findClass($ptr,$temp[0],\$foundHash);
		(UNIVERSAL::isa( $foundHash, "HASH" )) 
			or die "$0: Not found class $x at line $.\n";
		print STDERR "$#temp \n";
		if ($#temp>=1) {
			my $nameOfKind = $temp[1];
			
			my $kind = "variable";
			if ($nameOfKind=~s/\(\)$//) {
				$kind="function";
			}
			findKind($foundHash,$nameOfKind,\$foundHash,$kind);
			(UNIVERSAL::isa( $foundHash, "HASH" ))
				or die "$0: Not found function $x\n";
		}
		my $substitution = getDetailed($foundHash);
		s/\\ptexPaste\{([^\}]+)\}/$substitution/;
	}
	print;
}

sub getBrief
{
	my ($ptr)=@_;
	return $ptr->{"brief"};
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
		print "key=$item value=$value\n";
	}
}

sub findClass
{
	my ($ptr,$thisClass,$foundHash)=@_;
	for my $item (@$ptr) {
		if ($item->{"name"} eq $thisClass) {
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
#  		print "key=$item\n";

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
		}
		if ($item eq "name") {
			$lastName = $x;
		}
		
		if ($lastName eq $thisFunc and $lastKind eq $kind) {
			$$foundHash = $ptr;
# 			print STDERR "Found for hash $$foundHash\n";
			return;
		}
	}
}

sub findKindA
{
	my ($ptr,$thisFunc,$foundHash,$kind)=@_;
	for my $item (@$ptr) {
# 		print "key=$item\n";
		my $x = $item;
		if (UNIVERSAL::isa( $x, "HASH" )) {
			findKind($x,$thisFunc,$foundHash,$kind);
		} elsif (UNIVERSAL::isa( $x, "ARRAY" )) {
			findKindA($x,$thisFunc,$foundHash,$kind);
		} else {
# 			print "value=$x\n";
		}
	}
}

# sub procHash
# {
# 	my ($ptr)=@_;
# 	for my $item (keys %$ptr) {
# 		print "key=$item\n";
# 		my $x = $ptr->{$item};
# 		if (UNIVERSAL::isa( $x, "HASH" )) {
# 			procHash($x);
# 		} elsif (UNIVERSAL::isa( $x, "ARRAY" )) {
# 			procArray($x);
# 		} else {
# 			#print "value=$x\n";
# 		}
# 	}
# }
# 
# sub procArray
# {
# 	my ($ptr)=@_;
# 	for my $item (@$ptr) {
# 		#print "$item\n";
# 		my $x = $item;
# 		if (UNIVERSAL::isa( $x, "HASH" )) {
# 			procHash($x);
# 		} elsif (UNIVERSAL::isa( $x, "ARRAY" )) {
# 			procArray($x);
# 		} else {
# 			#print "value=$x\n";
# 		}
# 	}
# }
