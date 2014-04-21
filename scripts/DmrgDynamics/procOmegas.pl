#!/usr/bin/perl

use strict;
use warnings;
use lib "../ToolsDmrgDynamics";
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

sub createInput
{
	my ($n,$ind,$omega)=@_;
	my $file="input0.inp";
	open(FOUT,">$file") or die "$0: Cannot write to $file\n";
	my $steps = int($n/2) - 1;
	my $data = "data$ind.txt";
	my $nup = int($n/2);
	my $ndown = $nup;

	open(FILE,"$templateInput") or die "$0: Cannot open $templateInput: $!\n";

	while(<FILE>) {
		next if (/^#/);
		if (/\$([a-zA-Z0-9\[\]]+)/) {
				my $name = $1;
				my $str = "\$".$name;
				my $val = eval "$str";
				defined($val) or die "$0: Undefined substitution for $name\n";
				s/\$\Q$name/$val/g;
		}
		print FOUT;
	}

	close(FILE);
	close(FOUT);

	return $file;
}

sub createBatch
{
        my ($ind,$omega,$input) = @_;
        my $file = "Batch$ind.pbs";
        open(FOUT,">$file") or die "$0: Cannot write to $file: $!\n";

        open(FILE,"$templateBatch") or die "$0: Cannot open $templateBatch: $!\n";

        while(<FILE>) {
                while (/\$\$([a-zA-Z0-9\[\]]+)/) {
                        my $line = $_;
                        my $name = $1;
                        my $str = "\$".$name;
                        my $val = eval "$str";
                        defined($val) or die "$0: Undefined substitution for $name\n";
                        $line =~ s/\$\$$name/$val/;
                        $_ = $line;
                }

                print FOUT;
        }

        close(FILE);
        close(FOUT);

        print STDERR "$0: $file written\n";
        return $file;
}

sub submitBatch
{
        my ($batch) = @_;
        sleep(2);
        system("qsub $batch");
        print STDERR "$0: Submitted $batch\n";
}

sub execThis
{
	my ($cmd) = @_;
	print STDERR "$0: About to execute $cmd\n";
	system($cmd);
}

