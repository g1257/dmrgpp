#!/usr/bin/perl

use strict;
use warnings;
use Math::Trig;
use lib ".";
use OmegaUtils;

my ($templateInput,$templateBatch,$parallel) = @ARGV;
my $usage = "dollarizedInput dollarizedBatch howToSubmit\n";
$usage .="\t howToSubmit is one of nobatch  submit  test";
defined($parallel) or die "USAGE: $0 $usage\n";

my ($omega0,$total,$omegaStep,$obs,$GlobalNumberOfSites);
my $offset = 0;

my $isAinur = OmegaUtils::isAinur($templateInput);

my $hptr = {
"#OmegaBegin" => \$omega0,
"#OmegaTotal" => \$total,
"#OmegaStep" => \$omegaStep,
"#Observable" => \$obs,
"#Offset" => \$offset,
"TotalNumberOfSites" => \$GlobalNumberOfSites};

OmegaUtils::getLabels($hptr, $templateInput);

if ($omegaStep < 0) {
	my $beta = -$omegaStep;
	print STDERR "$0: Matsubara freq. assumed with beta= $beta\n";
	$omega0 = $omegaStep = 2.0*pi/$beta;
}

my $jobs = "";
my @outfiles;
for (my $i = $offset; $i < $total; ++$i) {
	my $omega = sprintf("%.3f", $omega0 + $omegaStep * $i);
	print STDERR "$0: About to run for omega = $omega\n";
	my ($jobid,$outfile) = runThisOmega($i, $omega, $obs, $isAinur);
	print STDERR "$0: Finished         omega = $omega\n";
	$jobs .= ":" if ($i > 0);
	$jobid =~ s/\..*$//;
	$jobs .= $jobid;
	push @outfiles, $outfile;
}

sub runThisOmega
{
	my ($ind, $omega, $obs, $isAinur) = @_;
	my $n = $GlobalNumberOfSites;
	my $input = createInput($n, $ind, $omega, $obs, $isAinur);
	my $jobid = "";
	my $outfile = "runFor$input";
	my $ext = ($isAinur) ? "ain" : "inp";
	$outfile =~ s/\.$ext//;
	$outfile .= ".cout";
	my $batch = createBatch($ind, $omega, $input, $obs);
	$jobid = submitBatch($batch, $parallel);
	#system("echo '#omega=$omega' >> $outfile") if ($submit eq "nobatch");
	return ($jobid,$outfile);
}

sub createInput
{
	my ($n, $ind, $omega, $obs, $isAinur)=@_;

	$n =~ s/;// if ($isAinur);
	my $ext = ($isAinur) ? "ain" : "inp";
	my $file="input$ind.$ext";
	open(FOUT, ">", "$file") or die "$0: Cannot write to $file\n";

	print FOUT "##Ainur1.0\n" if ($isAinur);

	my $steps = int($n/2) - 1;
	my $data = "data$ind.txt";
	my $nup = int($n/2);
	my $ndown = $nup;

	my %valuesHash = (
	"steps" => $steps,
	"data" => $data,
	"nup" => $nup,
	"ndown" => $ndown,
	"omega" => $omega,
	"obs" => $obs);

	open(FILE, "<", "$templateInput") or die "$0: Cannot open $templateInput: $!\n";

	while(<FILE>) {
		next if (/^#/);
		if (/\$([a-zA-Z0-9\[\]]+)/) {
				my $name = $1;
				my $val = $valuesHash{"$name"};
				defined($val) or die "$0: Undefined substitution for $name\n";
				s/\$\Q$name/$val/g;
		}

		print FOUT;
	}

	close(FILE);
	close(FOUT);

	return $file;
}

sub createFinalBatch
{
	my ($ind, $tarname, $files, $obs) = @_;
	createBatch($ind,0,"FINAL");
	my $fout = "temp.pbs";
	my $file = "Batch$ind.pbs";
	open(FOUT, ">", "$fout") or die "$0: Cannot write to $file: $!\n";
	open(FILE, "<", "$file") or die "$0: Cannot open $fout: $!\n";
	my @outfiles = @$files;
	while (<FILE>) {
		if (/FINAL/) {
			print FOUT "$tarname @outfiles\n";
			next;
		}

		print FOUT;
	}

	close(FOUT);
	close(FILE);
	system("cp $fout $file");
	unlink($fout);
	return $file;
}

sub createBatch
{
        my ($ind, $omega, $input, $obs) = @_;
        my $file = "Batch$ind.pbs";

	my %valuesHash = (
	"input" => $input,
	"ind" => $ind,
	"omega" => $omega,
	"obs" => $obs);

        open(FOUT, ">", "$file") or die "$0: Cannot write to $file: $!\n";

        open(FILE, "<", "$templateBatch") or die "$0: Cannot open $templateBatch: $!\n";

        while (<FILE>) {
                while (/\$\$([a-zA-Z0-9\[\]]+)/) {
                        my $line = $_;
                        my $name = $1;
                        my $val = $valuesHash{"$name"};
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
	my ($batch, $doIt) = @_;
	return if ($doIt ne "nobatch" and $doIt ne "submit");
	sleep(1);
	print STDERR "$0: Submitted $batch \n";
	my $execCommand = ($doIt eq "nobatch") ? "env PBS_O_WORKDIR=. /usr/bin/bash" : "qsub";
	my $ret = `$execCommand $batch`;
	chomp($ret);
	return $ret;
}


