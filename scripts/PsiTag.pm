#!/usr/bin/perl

use strict;
use warnings;
use utf8;

=pod
Syntax: tagging  mode content
content can appear multiline if scope is block
block is between ( ) multiline
tagging and mode must appear in the same line
automatic tagging is possible with .= content
Mode is either
             add definition and overwrite if it exits =
             add definition and append to it if it exits +=
             add definition or ignore if it exits =? 
             add definition or fail if it exits =!
             delete definition or ignore if it doesn't exist -=
             delete definition or fail if it exits -=!
 
tagging must not contain any of = + ? ! - < & *

Content is subject to interpretation as follows.
   
   (1) First non-whitespace character after mode in the same line, if it is a (

   (2) Last non-whitespace character in a line, if it is a )
 
   (3) First non-whitespace character in a line or after mode, if it is a <
       and what follows until newline is considered the tagging

   (4) First non-whitespace character in line or after mode, if it is an escape \
       will be removed. This is so that the next character will not be interpreted

   (5) Last non-whitespace character in a line, if it is a \
       will be removed. This is so that the previous character will not be interpreted

Tagging is canonicalized as follows: See canonicalTagName below.

=cut

package PsiTag;

sub main
{
	my ($tags, $file) = @_;

	my @lines;
	readLines(\@lines, $file);

	my %tags;
	readTags($tags, \@lines);
}

sub readLines
{
	my ($lines, $file) = @_;
	open(FILE, "<", $file) or die "$0: Cannot open $file : $!\n";
	while (<FILE>) {
		chomp;
		push @$lines, $_;
	}

	close(FILE);
}

sub readTags
{
	my ($tags, $lines) = @_;
	my $n = scalar(@$lines);
	my $blockScope = 0;
	my $multilineContent;
	my $multilineTag;
	my $multilineMode;
	for (my $i = 0; $i < $n; ++$i) {
		my $line = $lines->[$i];
		if (!$blockScope) {
			next if ($line =~ /^\#/ or isEmptyLine($line));
			if ($line =~ /^([^\=\+\?\!\-\<\&\*]+)([\=\+\?\!\-\<\&\*][^ \t]*)(.*)$/) {
				my $tag = $1;
				my $mode = $2;
				my $rest = $3;

				my $thisLineParens = 0;
				if ($rest =~ s/^[ \t]*\(//) { # (1)
					++$thisLineParens;
				}

				if ($rest =~ s/\)[ \t]*$//) { # (2) in line scope
					--$thisLineParens;
				}

				my $content = $rest."\n";
				$content =~ s/^[ \t]+//;
				

				if ($rest =~ /^[ \t]*\<(.*$)/) { # (3) in line scope
					my $existingTag = $1;
					$existingTag = canonicalTagName($existingTag);
					$content = $tags->{"$existingTag"}->{"content"};
					defined($content) or die "$0: No content for $existingTag\n";
				} else {
					$content =~ s/^[ \t]*\\//; # (4) in line scope
					$content =~ s/\\([ \t]*)$/$1/; # (5) in line scope
				}

				if ($thisLineParens < 0) {
					print STDERR "$0: ) found but not in block scope\n";
					syntaxError($line, $i + 1);
				}

				if ($thisLineParens > 0) {
					$blockScope = 1;
					$multilineContent = "";
					$multilineTag = canonicalTagName($tag);
					$multilineMode = $mode;
				} else {
					addTag($tags, canonicalTagName($tag), $mode, $content);
				}
			} else {
				syntaxError($line, $i + 1);
			}
		} else { # in block scope
			my $content = $line."\n";
			my $closeScope = 0;

			$content =~ s/^[ \t]+//;
			if ($content =~ s/\)[ \t]*$//) { # (2) in block scope
				$closeScope = 1;
			}

			if ($content =~ /^[ \t]*\<(.*$)/) { # (3) in block scope
				my $existingTag = $1;
				$existingTag = canonicalTagName($existingTag);
				my $ptr = $tags->{"$existingTag"};
				(ref($ptr) eq "HASH") or die "$0: Tag $existingTag not hash ref but ".ref($ptr)."\n";
				$content = $ptr->{"content"};
				defined($content) or die "$0: No content for $existingTag\n";
			} else {
				$content =~ s/^[ \t]*\\//; # (4) in block scope
				$content =~ s/\\([ \t]*)$/$1/; # (5) in block scope
			}

			$multilineContent .= $content;

			if ($closeScope) {
				$blockScope = 0;
				addTag($tags, $multilineTag, $multilineMode, $multilineContent);
				$multilineTag = $multilineMode = $multilineContent = "";
				next;
			}
		}
	}
}

sub addTag
{
	my ($tags, $tag, $mode, $content) = @_;
	#print STDERR "$0: Adding $tag, ignoring $mode\n";
	$tags->{"$tag"} = {"content" => $content};
}

sub syntaxError
{
	my ($line, $ind) = @_;
	print STDERR "$0: Syntax error line $ind\n";
	die "$0: ----> $line <------\n";
}

sub isEmptyLine
{
	my ($line) = @_;
	$line =~ s/[ \t]//g;
	return ($line eq "");
}

sub canonicalTagName
{
	my ($name) = @_;
	$name =~ s/\n/ /g;
	$name =~ s/^[ \t]*//;
	$name =~ s/[ \t]*$//;
	while ($name =~ s/[ \t][ \t]/ /) {}
	return $name;
}

1;


