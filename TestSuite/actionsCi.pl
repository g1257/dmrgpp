#!/usr/bin/perl

use strict;
use warnings;
use utf8;
use lib "../../scripts";
use timeObservablesInSitu;
use EnergyAncillaInSitu;
use CollectBrakets;
use Metts;

my ($action, $n, @what) = @ARGV;
defined($what[0]) or die "$0: USAGE: action n what\n";

my %actions = (getTimeObservablesInSitu => \&runTimeInSituObs,
	               getEnergyAncilla => \&runEnergyAncillaInSituObs,
	               CollectBrakets => \&runCollectBrakets,
	               metts => \&runMetts);

defined($actions{$action}) or die "$0: Action $action not registered\n";

$actions{$action}->($n, \@what);

sub runTimeInSituObs
{
	my ($n, $what) = @_;
	my $whatN = scalar(@$what);
	for (my $i = 0; $i < $whatN; ++$i) {
		my $file = "runForinput$n.cout";
		if (!(-r "$file")) {
			print STDERR "|$n|: WARNING: $file not readable\n";
			next;
		}

		my @temp = split(/ /, $what->[$i]);
		(scalar(@temp) == 2) or die "$0: FATAL annotation ".$what->[$i]."\n";
		my ($site, $label) = @temp;		
		my $fin;
		open($fin, "<", $file) or die "$0: Could not open $file : $!\n";
		my $fout;
		my $foutname = "timeObservablesInSitu${n}_$i.txt";
		if (!open($fout, ">", "$foutname")) {
			close($fin);
			die "$0: Could not write to $foutname: $!\n";
		}

		timeObservablesInSitu::main($site, $label, $fin, $fout);
		
		close($fin);
		close($fout);
	}
}

sub runEnergyAncillaInSituObs
{
	my ($n, $what) = @_;
	my $whatN = scalar(@$what);
	for (my $i = 0; $i < $whatN; ++$i) {
		my $file = "runForinput$n.cout";
		if (!(-r "$file")) {
			print STDERR "|$n|: WARNING: $file not readable\n";
			next;
		}

		my @temp = split(/ /, $what->[$i]);
		(scalar(@temp) == 2) or die "$0: FATAL annotation ".$what->[$i]."\n";
		my ($beta, $label) = @temp;
		my $fin;
		open($fin, "<", $file) or die "$0: Could not open $file : $!\n";
		my $fout;
		my $foutname = "energyAncillaInSitu${n}_$i.txt";
		if (!open($fout, ">", "$foutname")) {
			close($fin);
			die "$0: Could not write to $foutname: $!\n";
		}

		EnergyAncillaInSitu::main($beta, $label, $fin, $fout);

		close($fin);
		close($fout);
	}
}

sub runCollectBrakets
{
	my ($n, $what) = @_;
	my $whatN = scalar(@$what);
	for (my $i = 0; $i < $whatN; ++$i) {
		my $file = "runForinput$n.cout";
		if (!(-r "$file")) {
			print STDERR "|$n|: WARNING: $file not readable\n";
			next;
		}

		#my @temp = split(/ /, $what->[$i]); # arguments ot CollectBrakets
		my $foutname = "CollectBrakets${n}_$i.txt";

		CollectBrakets::main($file, $foutname);
	}
}

sub runMetts
{
	my ($n,$what) = @_;
	my $whatN = scalar(@$what);
	my %actions = ("Energy" => \&Metts::energy,
	               "Density" => \&Metts::density);
	for (my $i = 0; $i < $whatN; ++$i) {
		my $file = "runForinput$n.cout";
		if (!(-r "$file")) {
			print STDERR "|$n|: WARNING: $file not readable\n";
			next;
		}

		my @temp = split(/ /, $what->[$i]);
		(scalar(@temp) == 3) or next;

		my ($label, $arg0, $arg1) = @temp;
		if (($label ne "Energy") and ($label ne "Density")) {
			die "$0: Wrong annotation: $what->[$i]\n";
		}

		my $fin;
		open($fin, "<", $file) or die "$0: Could not open $file : $!\n";
		my $fout;
		my $foutname = "metts${n}_$i.txt";
		if (!open($fout, ">", "$foutname")) {
			close($fin);
			die "$0: Could not write to $foutname: $!\n";
		}

		my ($sum, $counter) = $actions{"$label"}($arg0, $arg1,0, $fin);
		print $fout "#"."$label=$sum $counter\n";

		close($fin);
		close($fout);
	}
}


