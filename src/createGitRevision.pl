#!/usr/bin/perl

use strict;
use warnings;

my ($file) = @ARGV;
defined($file) or exit(1);

my $hashPsimagLite = getGitHash("../../PsimagLite");

my $hashDmrgpp = getGitHash("..");

my $microArch = getMicroArch();

open(FOUT, ">", $file) or exit(2);
print FOUT "// DO NOT EDIT. It will be overwritten\n";
print FOUT "// Created by $0\n";
print FOUT "#define PSIMAGLITE_GIT_REV \"$hashPsimagLite\"\n";
print FOUT "#define DMRGPP_GIT_REV \"$hashDmrgpp\"\n";
print FOUT "#define MICRO_ARCH \"$microArch\"\n";
print FOUT "\n";
close(FOUT);

sub getGitHash
{
	my ($dir) = @_;
	my $file = "$dir/.git/HEAD";
	open(FILE, "<", $file) or return "E0";
	my $value = <FILE>;
	close(FILE);
	($value) or return "E1";
	chomp($value);
	if ($value =~ /^ref\: (.+$)/) {
		$file = "$dir/.git/$1";
	} else {
		return "E2";
	}

	open(FILE, "<", $file) or return "E3";
	$value = <FILE>;
	close(FILE);
	($value) or return "E4";
	chomp($value);

	my $code = system("cd $dir; git diff --quiet");
	($code) or $code = 0;
	$code >>= 8;
	$value .= " +M" if ($code != 0);

	return $value;
}

sub getMicroArch
{
	my $file = "/proc/cpuinfo";
	open(FILE, "<", $file) or return "E0";
	my $vendorId;
	while (<FILE>) {
		next unless (/vendor_id/);
		chomp;
		$vendorId = $_;
	}

	close(FILE);

	($vendorId) or return "E1";
	return "Intel" if ($vendorId =~ /intel/i);
	return "AMD" if ($vendorId =~ /AMD/);
	$vendorId =~ s/vendor_id[ \t] *:[ \t]*//;
	return $vendorId;
}


