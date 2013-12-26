#!/usr/bin/perl

use strict;
use warnings;

my ($file,$pbsname) = @ARGV;

defined($pbsname) or die "USAGE: $0 file pbsname\n";

my $ext=" ";

if ($file =~ /\.([^\.]*$)/) {
	$ext = $1;
}

if ($ext eq "inp") {
	makeTemplateInput($file);
} elsif ($ext eq "pbs") {
	makeTemplateBatch($file,$pbsname);
} else {
	die "$0: file must end in .inp or .pbs\n";
}

print STDERR "$0: WARNING: Check produced output before using\n";

sub makeTemplateBatch
{
	my ($file,$pbsname) = @_;

	defined($pbsname) or die "USAGE: $0 file.inp or $0 file.pbs pbsname\n";

	open(FILE,$file) or die "$0: Cannot open file $file : $!\n";
	while(<FILE>) {
		if (/#PBS +\-N +(.*$)/) {
			print "#PBS -N $pbsname\$\$site_\$\$site2\n";
			next;
		}

		print;
		last unless (/^#/);
	}

	close(FILE);

	print "cd \$PBS_O_WORKDIR\n";
	print "date\nperl dynamics\$\$dmrgOrLanczos.pl \$\$site \$\$site2 3 \$\$root \$\$operatorLabel\ndate\n\n";
}

sub makeTemplateInput
{
	my ($file) = @_;
	open(FILE,$file) or die "$0: Cannot open file $file : $!\n";

	while(<FILE>) {
		if (/^hubbardU/) {
			my @temp = split;
			die "$0: $_\n" if (scalar(@temp)<3);
			print "##U=$temp[2]\n";
			print "hubbardU \$hubbardU\n";
			next;
		}

		if (/^potentialV/) {
			my @temp = split;
			print "##V=$temp[2]\n";
			print "potentialV \$potentialV\n";
			next;
		}

		if (/^OutputFile=/) {
			print "OutputFile=\$data\n";
			next;
		}

		if (/^DynamicDmrgType=/) {
			print "DynamicDmrgType=\$type\n";
			next;
		}

		if (/^TSPSites/) {
			print "TSPSites \$sites\n";
			next;
		}

		if (/^TSPLoops/) {
			print "TSPLoops \$loops\n";
			next;
		}

		last if (/^TSPOperator=/);

		print;
	}

	close(FILE);

	print<<EOF;
TSPOperator=raw
RAW_MATRIX
\$hilbertSize \$hilbertSize
\$matrix[0]

TSPOperator=raw
RAW_MATRIX
\$hilbertSize \$hilbertSize
\$matrix[1]

TSPOperator=raw
RAW_MATRIX
\$hilbertSize \$hilbertSize
\$matrix[2]

TSPOperator=raw
RAW_MATRIX
\$hilbertSize \$hilbertSize
\$matrix[3]

IsPeriodicX=0
EOF

}


