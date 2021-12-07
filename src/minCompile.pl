#!/usr/bin/perl

use strict;
use warnings;
use utf8;
use lib ".";
use MinCompile;

my @models = @ARGV;

MinCompile::main("models", \@models);

