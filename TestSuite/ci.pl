#!/usr/bin/perl

use strict;
use warnings;
use Ci;

my ($min,$max,$submit) = @ARGV;
defined($submit) or $submit = "n";
my $templateBatch = "batchDollarized.pbs";
my @tests;
Ci::getTests(\@tests);

prepareDir();

my $total = scalar(@tests);
for (my $i = 0; $i < $total; ++$i) {
	my $n = $tests[$i];
	next if (defined($min) and $n < $min);
	next if (defined($max) and $n > $max);
	procTest($n,$submit);
}

sub procTest
{
	my ($n,$submit) = @_;
	my $cmd = "./dmrg -f ../inputs/input$n.inp";
	my $batch = createBatch($n,$cmd);
	submitBatch($batch) if ($submit eq "Y");
}

sub createBatch
{
	my ($ind,$cmd) = @_;
	my $file = "Batch$ind.pbs";
	open(FOUT,">$file") or die "$0: Cannot write to $file: $!\n";

	open(FILE,"../$templateBatch") or die "$0: Cannot open ../$templateBatch: $!\n";

	while (<FILE>) {
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
	my ($batch,$extra) = @_;
	defined($extra) or $extra = "";
	sleep(1);
	print STDERR "$0: Submitted $batch $extra $batch\n";

	my $ret = `qsub $extra $batch`;
	chomp($ret);
	return $ret;
}

sub prepareDir
{
	my $b = (-r "tests");
	system("mkdir tests") if (!$b);
	my $cmd = "diff ../src/dmrg tests/dmrg &> /dev/null";
	my $ret = system($cmd);
	system("cp -a ../src/dmrg tests/") if ($ret != 0);
	chdir("tests/");
}

