#!/usr/bin/perl -w

use strict;

my @uvector;
my ($vx,$vy);
my ($tx,$ty);
my $FileTst;
my $FileData;
my $BatchName;
my $InputFile;
my @ignoreVars=("PBS_O_WORKDIR");
	
my $tstSites=" 6 4 ";
my $templateInput="inputDollar.inp";
my $templateBatch="batchDollar.pbs";
my $rootDataFile="data";
my $rootTstFile="tst";
my $rootInput="input";
my $rootBatchName="dmrg14_";
my $rootBatchFile="batch";


my ($runs)=@ARGV;

if (defined($runs)) {
	askToSubmit($runs,"y");
	exit(0);
}

$vx=$vy=0.0;
my $n = readLabel($templateInput,"TotalNumberOfSites=");

fillVector(\@uvector,$n,10);

my $counter=0;
for ($tx=0.0;$tx<1.1;$tx+=0.1) {
	$ty = 1.0-$tx;
	$ty = 0 if ($ty<1e-6);
	createRun($counter);
	createBatch($counter);
	$counter++;
}

askToSubmit($counter,"n");

sub createRun
{
	my ($counter)=@_;
	
	$FileTst="$rootTstFile$counter.txt";
	$FileData="$rootDataFile$counter.txt";
	
	
	my $fout = "$rootInput$counter.inp";
	unDollarize($templateInput,$fout);
}

sub createBatch
{
	my ($counter)=@_;
	
	$BatchName="$rootBatchName$counter";
	$InputFile="$rootInput$counter.inp";
	my $batchFile="$rootBatchFile$counter.pbs";
	
	unDollarize($templateBatch,$batchFile);
}

sub readLabel
{
	my ($file,$label)=@_;
	my $ret;
	open(FILE,$file) or die
		"Cannot open $file: $!\n";
	while(<FILE>) {
		chomp;
		if (/^$label(.*$)/) {
			$ret=$1;
			last;
		}
	}
	close(FILE);
	defined($ret) or die "readLabel: Not found $label in $file\n";
	return $ret;
}

sub fillVector
{
	my ($v,$n,$val)=@_;
	for (my $i=0;$i<$n;$i++) {
		$v->[$i] = $val;
	}
}



sub unDollarize
{
	my ($inFile,$outFile)=@_;
	open(FILE,$inFile) or die
		"Cannot open $inFile: $!\n";
	open(FOUT,">$outFile") or die
		"Cannot open $outFile for writing: $!\n"; 
	while(<FILE>) {
		my $line = $_;
		while (s/([\$\@])([a-zA-Z_]+)//) {
			my $indicator = $1;
			my $var=$2;
			next if (isInVector(\@ignoreVars,$var));
			
			my @tmp=eval("$indicator$var");
			defined($tmp[0]) or die "Undefined for $var\n";
			#print STDERR "$indicator$var @tmp\n";
			if ($indicator eq "\$") {
				$line =~ s/\$$var/@{tmp}/g;
			} else {
				$line =~ s/\@$var/@{tmp}/g;
			}
		}
		print FOUT "$line";
	}
	close(FILE);
	close(FOUT);
}

sub isInVector
{
	my ($array,$var)=@_;
	foreach my $tmp (@$array) {
		return 1 if ($var eq $tmp);
	}
	return 0;
}

sub askToSubmit
{
	my ($total,$d)=@_;
	print "Do you want to submit all jobs now?\n";
	print "Available: y or n\n";
	print "Default is: $d (press ENTER)";
	$_=<STDIN>;
	chomp;
	$_="$d" if ($_ eq "");
	return if ($_=~/n/i);
	submitAll($total);
	
}

sub submitAll
{
	my ($total)=@_;
	for (my $i=0;$i<$total;$i++) {
		my $batchFile="$rootBatchFile$i.pbs";
		system("qsub $batchFile");
	}
}

