#!/usr/bin/perl

use strict;
use warnings;
use utf8;

my $makeJ = shift @ARGV;
defined($makeJ) or die "USAGE: $0 makeJ\n";

main(4, $makeJ, \@ARGV);

sub main
{
	my ($items, $makeJ, $codes) = @_;

	for (my $i = 0; $i < $items; ++$i) {
		kompileRig($i, $makeJ, $codes);
	}

	print "---------------------\n";
}

sub kompileRig
{
	my ($ind, $makeJ, $codes) = @_;
	my $command = "make clean; make -j $makeJ";
	my @paths = qw!../../PsimagLite/lib !;
	addCodes(\@paths, $codes);
	my $n = scalar(@paths);
	for (my $i = 0; $i < $n; ++$i) {
		kompileRigEach($ind, $paths[$i], $command, $codes->[$i]);
	}
}

sub kompileRigEach
{
	my ($ind, $path, $command, $code) = @_;
	my $psiTag = "../../dmrgpp/TestSuite/inputs/KompileRig.psiTag";
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
	my ($paths, $codes) = @_;
	addAll($codes);
	my $n = scalar(@$codes);
	for (my $i = 0; $i < $n; ++$i) {
		my $code = $codes->[$i];
		if ($code eq "dmrgpp") {
			push @$paths, "../src";
		} elsif ($code eq "LanczosPlusPlus") {
			push @$paths, "../../LanczosPlusPlus/src";
		} elsif ($code eq "BetheAnsatz") {
			push @$paths, "../../BetheAnsatz/src";
		} elsif ($code eq "FreeFermions") {
			push @$paths, "../../FreeFermions/examples";
		} elsif ($code eq "merapp") {
			push @$paths, "../../merapp/src";
		} elsif ($code eq "PsimagLite") {
			push @$paths, " ../../PsimagLite/drivers";
		} else {
			die "$0: Unknown code $code\n";
		}
	}
}

sub addAll
{
	my ($codes) = @_;
	return if (scalar(@$codes) > 0);
	@$codes = qw/PsimagLite dmrgpp  LanczosPlusPlus/; # BetheAnsatz  FreeFermions merapp/;
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

