#!/usr/bin/perl

use strict;
use warnings;
use utf8;

=pod
Syntax: tagging  mode content
content can appear multiline if scope is block
block is between ( ) multiline
tagging and mode must appear in the same line

Mode is either
             add definition and overwrite if it exits =
             add definition and append to it if it exits +=
             add definition or ignore if it exits =?
             add definition or fail if it exits =!
             delete from definition or ignore if it doesn't exist -=
             delete from definition or fail if it exits -=!

tagging must not contain any of = + ? ! - < & *

Content is subject to interpretation as follows.

   (0) sometag = ()
	Deletes sometag if it exits. Spaces around =  are optional. Note that
       sometag =
       sets the content of sometag to nothing, and is different from deletion.

   (1) First non-whitespace character after mode in the same line, if it is a (

   (2) Last non-whitespace character in a line, if it is a )

   (3) First non-whitespace character in a line or after mode, if it is a <
       and what follows until newline is considered the tagging

   (4) First non-whitespace character in line or after mode, if it is an escape \
       will be removed. This is so that the next character will not be interpreted

   (5) Last non-whitespace character in a line, if it is a \
       will be removed. This is so that the previous character will not be interpreted

Tagging is canonicalized as follows: See canonicalTagName below.

Proposal for improvement: if statement after mode

Example:

mytag ? [if xyz] = 42

sets mytag to 42 only if xyz is an existing tag

Syntax: tagging ? [ modifier ] mode content

where tagging, mode, and content are treated and defined as before,
modifier MUST NOT be multiline; and \[ is replaced into [ and \] into ]

modifer = submodifer1 &&& submodifer2 &&& submodifer3

\&&& is replaced into &&& for each submodifer.

The only submodifer type allowed for now is type if. And it is allowed only once.

submodifierOfTypeIf = if predicate

where

predicate = subpredicate1 SEPARATOR subpredicate2 SEPARATOR ...

If only one subpredicate is present SEPARATOR MUST be ommited.
Separator is one of || for logical OR, and && for logical AND.
For each subpredicate, \|| is replaced into || and \&& into &&.

subpredicate = tag
subpredicate = ! tag
subpredicate = tag == content
subpredicate = tag != content
subpredicate = tag1 == < tag2
subpredicate = tag1 != < tag2
Spaces are optional. Spaces before and after content will be part of content.
content MUST NOT be multiline; the whole predicate is within the line.

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
	my $parensBalance = 0;

	for (my $i = 0; $i < $n; ++$i) {
		my $line = $lines->[$i];
		if ($line =~ /^[ \t]*\}[ \t]*$/) {
			print STDERR "$0: Warning Closing brace on its own in line $i\n";
		}

		if (!$blockScope) {

			next if ($line =~ /^\#/ or isEmptyLine($line));

			if ($line =~ /^([^\=\+\?\!\-\<\&\*]+)\=[ \t]*\(\)[ \t]*$/) {
				my $tag = $1;
				removeTag($tags, canonicalTagName($tag));
				next;
			}

			if ($line =~ /^([^\=\+\?\!\-\<\&\*]+)([\=\+\?\!\-\<\&\*]+)[ \t]*(.*)$/) {
				my $tag = $1;
				my $mode = $2;
				my $rest = $3;

				my $thisLineParens = 0;
				if ($rest =~ s/^[ \t]*\(//) { # (1)
					++$thisLineParens;
					++$parensBalance;
					die "$0: FATAL: (line scope) Nested parens not allowed, line=$line\n" if ($parensBalance > 1);
				}

				if ($rest =~ s/\)[ \t]*$//) { # (2) in line scope
					--$thisLineParens;
					--$parensBalance;
					die "$0: FATAL: Closing parens but context closed already, line=$line\n" if ($parensBalance < 0);
				}

				my $content = $rest."\n";
				$content =~ s/^[ \t]+//;

				if ($thisLineParens < 0) {
					print STDERR "$0: ) found but not in block scope\n";
					syntaxError($line, $i + 1);
				}

				($tag) or die "$0: FATAL: tag does not exist in line scope $line\n";
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
				--$parensBalance;
				die "$0: FATAL: (block scope) Nested parens not allowed, line=$line\n" if ($parensBalance > 1);
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
	# mode is 
	# = (add or overwrite) 
	# += (append) 
	# -= (delete) 
	# =? (selective add)
	# =! (strict add)
	# -=! (strict delete)
	# Does the tag exist?
	my $mytag = $tags->{"$tag"};
	my $b = defined($mytag);
	my $oldContent = ($b) ? $tags->{"$tag"}{"content"} : "";
	my $append = (($mode eq "=?" && !$b) or ($mode eq "+=") or ($mode eq "=!" && !$b) or ($mode eq "="));
	
	die "$0: FATAL: Adding tag $tag with $mode, but $tag already exists\n" if ($mode eq "=!" && $b);

	return if ($mode eq "=?" && $b);

	if ($append) {
		$oldContent = "" if ($mode eq "=");
		$tags->{"$tag"} = {"content" => "$oldContent$content"};
		return;
	}

	if ($mode eq "-=" or $mode eq "-=!") {

		print STDERR "$0: WARNING: Deletion of content in non-existing tag $tag\n" if (!$b);

		my $hasIt = ($oldContent =~ /\Q$content/);

		die "$0: FATAL: Strict delete of non-existing content in tag $tag\n" if ($mode eq "-=!" && !$hasIt);

		return if (!$hasIt);

		$oldContent =~ s/\Q$content//g if ($hasIt);
		$tags->{"$tag"} = {"content" => "$oldContent"};
		return;
	}

	die "$0: FATAL: Unknown mode $mode applied to tag $tag\n";	
}

sub removeTag
{
	my ($tags, $tag) = @_;

	defined($tags->{"$tag"}) or return;	
	delete $tags->{"$tag"};
}

sub syntaxError
{
	my ($line, $ind) = @_;
	print STDERR "$0: FATAL: Syntax error line $ind\n";
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

sub unWrap
{
	my ($tags, $text) = @_;
	my @lines = split/\n/, $text;
	my $n = scalar(@lines);
	my $result = "";

	for (my $i = 0; $i < $n; ++$i) {
		my $line = $lines[$i];
		my $content = $line;
		if ($line =~ /^[ \t]*\<(.*$)/) { # (3) in block scope
			my $existingTag = $1;
			($existingTag) or die "$0: FATAL: Tag does not exist in unWrap: $line\n";
			$existingTag = canonicalTagName($existingTag);
			my $ptr = $tags->{"$existingTag"};
			defined($ptr) or die "$0: FATAL: Tag $existingTag doesn't exist\n";
			(ref($ptr) eq "HASH") or die "$0: FATAL: Tag $existingTag not hash ref but ".ref($ptr)."\n";
			$content = $ptr->{"content"};
			defined($content) or die "$0: FATAL: No content for $existingTag\n";
			$content = unWrap($tags, $content);
		}

		$content =~ s/^[ \t]*\\//; # (4)
		$content =~ s/\\([ \t]*)$/$1/; # (5)

		$result .= $content."\n";
	}

	return $result;
}

1;


