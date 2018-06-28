#!/usr/bin/perl

use strict;
use warnings;
use utf8;

use lib ".";
use Combinatorial;

my @levels;
my $counter = 0;
$levels[$counter++] = ["", "CPPFLAGS += -DUSE_FLOAT"];
$levels[$counter++] = ["", "CXX = clang++ -mtune=native"];

my ($makeJ) = @ARGV;
defined($makeJ) or die "USAGE: $0 makeJ\n";

Combinatorial::main(\@levels, \&myCallback);

sub myCallback
{
	my ($items) = @_;
	
	print flattenWithNewLines($items);

	kompileRig($items);

	print "---------------------\n";
}

sub kompileRig
{
	my ($items) = @_;
	my $textToAppend = flattenWithNewLines($items);
	my $command = "make clean; make -j $makeJ";
	my @paths = ("../../PsimagLite/lib", "../../PsimagLite/drivers", "../src");
	my $n = scalar(@paths);
	for (my $i = 0; $i < $n; ++$i) {
		kompileRigEach($paths[$i], $command, $textToAppend);
	}
}

sub kompileRigEach
{
	my ($path, $command, $textToAppend) = @_;
	unlink($path."/Config.make");
	my $cmd = "cd $path; echo Y | perl configure.pl";
	executeAndDieIfNotSuccess($cmd);
	appendToFile($path."/Config.make", $textToAppend);
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

