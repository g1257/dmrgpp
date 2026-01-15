#!/usr/bin/perl
=pod
Copyright (c) 2025 UT-Battelle, LLC
All rights reserved

PD I provide this file to allow the cmake to call the
DmrgDriver::createTemplates
in the simplest way possible to avoid run_drivers.pl.

*********************************************************
THE SOFTWARE IS SUPPLIED BY THE COPYRIGHT HOLDERS AND
CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED
WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A
PARTICULAR PURPOSE ARE DISCLAIMED.

Please see full open source license included in file LICENSE.
*********************************************************

=cut
use warnings;
use strict;

use DmrgDriver;

# This returns an array of 0's length number of files that would be
# created.
my $templates = DmrgDriver::createTemplates(1);

# compare $i to the number of members of templates
