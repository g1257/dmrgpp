#!/usr/bin/perl

use strict;
use warnings;
use utf8;

package PsiTag;

sub main
{
	my ($file) = @_;

	my @lines;
	readLines(\@lines, $file);

	my %tags;
	readTags(\%tags, \@lines);
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
				if ($rest =~ s/^[ \t]*\(//) {
					++$thisLineParens;
				}

				if ($rest =~ s/\)[ \t]*$//) {
					--$thisLineParens;
				}

				my $content = $rest;
				if ($rest =~ /^[ \t]*\<(.*$)/) {
					my $existingTag = $1;
					$existingTag = canonicalTagName($existingTag);
					$content = $tags->{"$existingTag"}->{"content"};
					defined($content) or die "$0: No content for $existingTag\n";
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
			my $content = $line;
			my $closeScope = 0;

			if ($content =~ s/\)[ \t]*$//) {
				$closeScope = 1;
			}

			if ($content =~ /^[ \t]*\<(.*$)/) {
				my $existingTag = $1;
				$existingTag = canonicalTagName($existingTag);
				my $ptr = $tags->{"$existingTag"};
				(ref($ptr) eq "HASH") or die "$0: Tag $existingTag not hash ref but ".ref($ptr)."\n";
				$content = $ptr->{"content"};
				defined($content) or die "$0: No content for $existingTag\n";
			}

			if ($closeScope) {
				$blockScope = 0;
				addTag($tags, $multilineTag, $multilineMode, $multilineContent);
				$multilineTag = $multilineMode = $multilineContent = "";
				next;
			}

			$multilineContent .= $content;
		}
	}
}

sub addTag
{
	my ($tags, $tag, $mode, $content) = @_;
	print STDERR "$0: Adding $tag, ignoring $mode\n";
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
	$name =~ s/^[ \t]*//;
	$name =~ s/[ \t]*$//;
	while ($name =~ s/[ \t][ \t]/ /) {}
	return $name;
}

1;


