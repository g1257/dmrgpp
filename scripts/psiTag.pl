#!/usr/bin/perl

use strict;
use warnings;
use utf8;
use lib ".";
use PsiTag;

my ($file) = @ARGV;

PsiTag::main($file);

