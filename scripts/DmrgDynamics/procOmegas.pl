#!/usr/bin/perl

use strict;
use warnings;
use Utils;

my ($omega0,$total,$omegaStep,$what,$q) = @ARGV;
$q = 0 if (defined($what) and $what eq "optical");
defined($q) or die "USAGE: $0 omegaBegin omegaTotal omegaStep what q\n";

my $templateInput = "inputTemplate.inp";
my $templateBatch = "batchTemplate.pbs";

for (my $i = 0; $i < $total; ++$i) {
	my $omega = $omega0 + $omegaStep * $i;
	print STDERR "$0: About to proc for omega = $omega\n";
	procThisOmega($i,$omega);
	print STDERR "$0: Finished         omega = $omega\n";
}

print STDERR "$0: Final result in out.txt\n";

sub procThisOmega
{
	my ($ind,$omega) = @_;
	my $n = Utils::getLabel($templateInput,"TotalNumberOfSites=");
	my $cmd = "perl correctionVector.pl $omega < out$ind.txt > out$ind.space";
	execThis($cmd);
	$cmd = "perl ft.pl < out$ind.space > out$ind.sq";
	execThis($cmd) if ($what eq "freq");

	if ($what eq "space") {
		my $site = $q;
		$cmd = "perl extractQ.pl $site  out$ind.space >> out.txt";
	} elsif ($what eq "freq") {
		$cmd = "perl extractQ.pl $q  out$ind.sq >> out.txt";
	} elsif ($what eq "optical") {
		$cmd = "perl optical.pl < out$ind.space >> out.txt";
	} else {
		die "$0: Unknown action $what\n";
	}

	execThis($cmd);
}

sub execThis
{
	my ($cmd) = @_;
	print STDERR "$0: About to execute $cmd\n";
	system($cmd);
}

