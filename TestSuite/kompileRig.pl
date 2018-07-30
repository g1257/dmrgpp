#!/usr/bin/perl

use strict;
use warnings;
use utf8;

my ($makeJ) = @ARGV;
defined($makeJ) or die "USAGE: $0 makeJ\n";

main(4, $makeJ);

sub main
{
	my ($items, $makeJ) = @_;

	for (my $i = 0; $i < $items; ++$i) {
		kompileRig($i, $makeJ);
	}

	print "---------------------\n";
}

sub kompileRig
{
	my ($ind, $makeJ) = @_;
	my $command = "make clean; make -j $makeJ";
	my @paths = ("../../PsimagLite/lib", "../../PsimagLite/drivers", "../src");
	my $n = scalar(@paths);
	for (my $i = 0; $i < $n; ++$i) {
		kompileRigEach($ind, $paths[$i], $command);
	}
}

sub kompileRigEach
{
	my ($ind, $path, $command) = @_;
	my $psiTag = "../../dmrgpp/TestSuite/inputs/KompileRig.psiTag";
	my $cmd = "cd $path; perl configure.pl -f KompileRig$ind -c $psiTag";
	executeAndDieIfNotSuccess($cmd);
	$cmd = "cd $path; $command";
	executeAndDieIfNotSuccess($cmd);
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

