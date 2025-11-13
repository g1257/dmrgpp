#!/usr/bin/perl

use strict;
use warnings;
use utf8;
use lib ".";
use PsiTag;

my ($file, $label) = @ARGV;
die "USAGE: $0 file [label]\n" if (!defined($file));

my %tags;
PsiTag::main(\%tags, $file);

exit(0) if (!defined($label));

my $ptr = $tags{"$label"};

defined($ptr) or die "$0: No tag named $label\n";

my $content = $ptr->{"content"};

defined($content) or die "$0: No content for tag $label\n";

print $content;


