#!/usr/bin/perl

use strict;
use warnings;

use lib '../scripts';
use Make;

my @drivers = qw/integrator sparseSolverTest testCRSMatrix minimizer
                 combineContinuedFraction continuedFractionCollection range
                 kernelPolynomial linearPrediction options randomTest svd
		 testLapack threads testIsClass testMemResolv1 sumDecomposition
		 calculator closuresTest base64test testIoNg/;
my $fh;
open($fh, ">", "Makefile") or die "$0: Cannot open Makefile for writing\n";

Make::newMake($fh, \@drivers, {"code" => "PsimagLite"});

close($fh);

print "$0: Makefile has been written\n";

