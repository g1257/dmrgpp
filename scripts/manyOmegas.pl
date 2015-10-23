#!/usr/bin/perl

use strict;
use warnings;
use Math::Trig;
use OmegaUtils;

my ($templateInput,$templateBatch,$parallel) = @ARGV;
my $usage = "dollarizedInput dollarizedBatch howToSubmit\n";
$usage .="\t howToSubmit is one of nobatch  submit  test";
defined($parallel) or die "USAGE: $0 $usage\n";

my ($omega0,$total,$omegaStep,$obs,$GlobalNumberOfSites);
my $offset = 0;

my $hptr = {
"#OmegaBegin" => \$omega0,
"#OmegaTotal" => \$total,
"#OmegaStep" => \$omegaStep,
"#Observable" => \$obs,
"#Offset" => \$offset,
"TotalNumberOfSites" => \$GlobalNumberOfSites};

OmegaUtils::getLabels($hptr,$templateInput);

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
	my ($jobid,$outfile) = runThisOmega($i,$omega);
	print STDERR "$0: Finished         omega = $omega\n";
	$jobs .= ":" if ($i > 0);
	$jobid =~ s/\..*$//;
	$jobs .= $jobid;
	push @outfiles, $outfile;
}

my $tarname = "runFor".$templateInput;
$tarname =~ s/\.inp//;
$tarname .= ".tar.gz";
$tarname = "tar --group nobody --owner nobody -zcvf $tarname";
if ($parallel eq "nobatch") {
	system("$tarname @outfiles");
} else {
	my $batch = createFinalBatch($total,$tarname,\@outfiles);
	my $afterok =  "-W depend=afterok:$jobs";
	submitBatch($batch,$afterok) if ($parallel eq "submit");
}

sub runThisOmega
{
	my ($ind,$omega) = @_;
	my $n = $GlobalNumberOfSites;
	my $input = createInput($n,$ind,$omega);
	my $jobid = "";
	my $outfile = "runFor$input";
	$outfile =~ s/\.inp//;
	$outfile .= ".cout";
	if ($parallel eq "nobatch") {
		(-x "./dmrg") or die "$0: No ./dmrg found in this directory\n";
		system("./dmrg -f $input -o ':$obs.txt'");
		system("echo '#omega=$omega' >> $outfile");
	} else {
		my $batch = createBatch($ind,$omega,$input);
		$jobid = submitBatch($batch) if ($parallel eq "submit");
	}

	return ($jobid,$outfile);
}

sub createInput
{
	my ($n,$ind,$omega)=@_;
	my $file="input$ind.inp";
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

sub createFinalBatch
{
	my ($ind,$tarname,$files) = @_;
	createBatch($ind,0,"FINAL");
	my $fout = "temp.pbs";
	my $file = "Batch$ind.pbs";
	open(FOUT,">$fout") or die "$0: Cannot write to $file: $!\n";
	open(FILE,"$file") or die "$0: Cannot open $fout: $!\n";
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
        my ($ind,$omega,$input) = @_;
        my $file = "Batch$ind.pbs";
        open(FOUT,">$file") or die "$0: Cannot write to $file: $!\n";

        open(FILE,"$templateBatch") or die "$0: Cannot open $templateBatch: $!\n";

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

