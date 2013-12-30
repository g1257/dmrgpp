#!/usr/bin/perl

use strict;
use warnings;
use Utils;

my ($dmrgOrLanczos,$root,$operatorLabel,$useReflectionSymmetry,$submit) = @ARGV;
defined($useReflectionSymmetry) or die "USAGE: $0 dmrgOrLanczos root operatorLabel useReflectionSymmetry submit\n";
Utils::checkRange($dmrgOrLanczos,"Lanczos","Dmrg");
Utils::checkRange($operatorLabel,"c","n","s+","s-","sz");
Utils::checkRange($submit,"nobatch","nosubmit","submit");

my $templateInput = "inputTemplate.inp";
my $templateBatch = "batchTemplate.pbs";
my $n = Utils::getLabel($templateInput,"TotalNumberOfSites=");

for (my $site=0; $site<$n; $site++) {
	for (my $site2=$site; $site2<$n; $site2++) {
		next if ($useReflectionSymmetry and Utils::reflected($site,$site2,$n));
		
		if ($submit eq "nobatch") {
			noBatch($site,$site2,$dmrgOrLanczos,$root,$operatorLabel);
			next;
		}

		my $batch = createBatch($site,$site2,$dmrgOrLanczos,$root,$operatorLabel);
		system("sync");
		submitBatch($batch) if ($submit eq "submit");
	}
}

sub noBatch
{
	my ($site,$site2,$dmrgOrLanczos,$root,$operatorLabel) = @_;
	system("perl dynamics$dmrgOrLanczos.pl $site $site2 3 $root $operatorLabel");
}


sub createBatch
{
	my ($site,$site2,$dmrgOrLanczos,$root,$operatorLabel) = @_;
	my $ind = "$site"."_"."$site2";
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


