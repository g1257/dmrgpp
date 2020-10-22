#!/usr/bin/perl

use strict;
use warnings;
use utf8;

my $makeJ = shift @ARGV;
defined($makeJ) or die "USAGE: $0 makeJ\n";

main(4, $makeJ, \@ARGV);

sub main
{
	my ($items, $makeJ, $clineCodes) = @_;

	for (my $i = 0; $i < $items; ++$i) {
		kompileRig($i, $makeJ, $clineCodes);
	}

	print "---------------------\n";
}

sub kompileRig
{
	my ($ind, $makeJ, $clineCodes) = @_;
	my $command = "make clean; make -j $makeJ";
	my @paths;
	my @codes = addCodes(\@paths, $clineCodes);
	my $n = scalar(@paths);
	die "$0: Number of Codes and Paths must be equal but $n != ".scalar(@codes)."\n"
	if ($n != scalar(@codes));
	for (my $i = 0; $i < $n; ++$i) {
		kompileRigEach($ind, $paths[$i], $command, $codes[$i]);
	}
}

sub kompileRigEach
{
	my ($ind, $path, $command, $code) = @_;
	return if ($code eq "cincuenta" and $ind < 2);

	print STDERR "-------------->> Compiling $code with profile $ind\n";
	my $psiTag0 = "TestSuite/inputs/KompileRig.psiTag";
	my $psiTag = "../../dmrgpp/$psiTag0";
	$psiTag = $psiTag0 if (-r "$path/$psiTag0");
	my $cmd = "cd $path; perl configure.pl -f KompileRig$ind -c $psiTag";
	executeAndDieIfNotSuccess($cmd);
	$cmd = "cd $path; $command";
	executeAndDieIfNotSuccess($cmd);
	return if ($code ne "dmrgpp");
	$cmd = "cd $path; cd ../doc; make manual.pdf";
	executeAndDieIfNotSuccess($cmd);
}

sub addCodes
{
	my ($paths, $clineCodes) = @_;
	my @codes = addAll($paths, $clineCodes);
	my $n = scalar(@codes);
	for (my $i = 0; $i < $n; ++$i) {
		my $code = $codes[$i];
		if ($code eq "PsimagLiteLib") {
		} elsif ($code eq "dmrgpp") {
			push @$paths, "../src";
		} elsif ($code eq "LanczosPlusPlus") {
			push @$paths, "../../LanczosPlusPlus/src";
		} elsif ($code eq "BetheAnsatz") {
			push @$paths, "../../BetheAnsatz/src";
		} elsif ($code eq "FreeFermions") {
			push @$paths, "../../FreeFermions/examples";
		} elsif ($code eq "merapp") {
			push @$paths, "../../merapp/src";
		} elsif ($code eq "PsimagLiteDrivers") {
			push @$paths, " ../../PsimagLite/drivers";
		} elsif ($code eq "cincuenta") {
			push @$paths, "../../cincuenta/src";
		} else	{
			die "$0: Unknown code $code\n";
		}
	}

	return @codes;
}

sub addAll
{
	my ($paths, $clineCodes) = @_;
	push @$paths, "../../PsimagLite/lib ";
	my @codes = ("PsimagLiteLib");
	my $n = scalar(@$clineCodes);
	if ($n > 0) {
		push @codes, @$clineCodes;
		return @codes;
	}

	my @codes2 = qw/PsimagLiteDrivers dmrgpp  LanczosPlusPlus cincuenta/; # BetheAnsatz  FreeFermions merapp/;
	push @codes, @codes2;
	return @codes;
}

sub flattenWithNewLines
{
	my ($items) = @_;
	my $text = "";
	my $n = scalar(@$items);
	for (my $i = 0; $i < $n; ++$i) {
		$text .= "$items->[$i]\n";
	}

	return $text;
}

sub appendToFile
{
	my ($file, $textToAppend) = @_;
	open(FILE, ">>", $file) or die "$0: Cannot open $file : $!\n";
	print FILE $textToAppend;
	close(FILE);
}

sub executeAndDieIfNotSuccess
{
	my ($cmd) = @_;
	system($cmd);
	die "$0: Command $cmd FAILED\n" if ($? != 0);
}

