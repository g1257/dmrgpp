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
my $hptr = {"#OmegaBegin" => \$omega0,
            "#OmegaTotal" => \$total,
			"#OmegaStep" => \$omegaStep,
			"#Observable" => \$obs,
			"TotalNumberOfSites" => \$GlobalNumberOfSites};

OmegaUtils::getLabels($hptr,$templateInput);

if ($omegaStep < 0) {
	my $beta = -$omegaStep;
	print STDERR "$0: Matsubara freq. assumed with beta= $beta\n";
	$omega0 = $omegaStep = 2.0*pi/$beta;
}

for (my $i = 0; $i < $total; ++$i) {
	my $omega = sprintf("%.3f", $omega0 + $omegaStep * $i);
	print STDERR "$0: About to run for omega = $omega\n";
	runThisOmega($i,$omega);
	print STDERR "$0: Finished         omega = $omega\n";
}

sub runThisOmega
{
	my ($ind,$omega) = @_;
	my $n = $GlobalNumberOfSites;
	my $input = createInput($n,$ind,$omega);
	if ($parallel eq "nobatch") {
		(-x "./dmrg") or die "$0: No ./dmrg found in this directory\n";
		system("./dmrg -f $input -o ':$obs.txt' &> out$ind.txt");
		system("echo '#omega=$omega' >> out$ind.txt");
	} else {
		my $batch = createBatch($ind,$omega,$input);
		submitBatch($batch) if ($parallel eq "submit");
	}
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
        my ($batch) = @_;
        sleep(2);
        system("qsub $batch");
        print STDERR "$0: Submitted $batch\n";
}

