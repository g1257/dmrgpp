#!/usr/bin/perl
#
use strict;
use warnings;
use utf8;

my ($uroot, $lOrUp) = @ARGV;

defined($uroot) or die "USAGE: $0 root [L | U | sz]\n";

if (!defined($lOrUp)) {
	my $file = "sample.tex";
	my $fout = "$uroot.tex";
	my $name = "$uroot.pgfplots";
	fromTexToTex($fout, $file, $name);;
	exit(0);
}

my $root = "$uroot$lOrUp"."_ky";

doFile(0, $root);
doFile(1, $root);

sub doFile
{
	my ($ind, $root) = @_;

	my $file = "outSpectrum$ind.pgfplots";
	return unless (-r "$file");
	my $fout = "$root$ind.pgfplots";
	unlink("$fout");
	my $cmd = "cp $file $fout";
	die "$0: Failed to create $fout\n" if (-r "$fout");
	system("$cmd");

	$file = "sample.tex";
	$fout = "$root$ind.tex";

	my $name = $fout;
	$name =~ s/tex$/pgfplots/;
	fromTexToTex($fout, $file, $name);
}

sub fromTexToTex
{
	my ($fout, $file, $name) = @_;

	my $dirForTex = $0;
	$dirForTex =~ s/pgfplot\.pl$//;
	die "$0: Not a directory $dirForTex\n" unless (-d "$dirForTex");

	system("cp $dirForTex/sample.tex .");
	system("cp $dirForTex/palette.tex .");

	copyAndEdit($fout, $file, $name);

	my $cmd = "pdflatex $fout";
	system("$cmd");
	sleep(1);
	my $foutpdf = $fout;
	$foutpdf =~ s/\.tex$/.pdf/;
	$cmd = "/usr/bin/pdftocairo -singlefile -png $foutpdf";
	system("$cmd") if (-x "/usr/bin/pdftocairo");
}

sub copyAndEdit
{
	my ($fout, $file, $name) = @_;
	open(FILE, "<", "$file") or die "$0: Cannot open $file : $!\n";
	my $ret = open(FOUT, ">", "$fout");
	if (!$ret) {
		close(FILE);
		die "$0: Cannot open write to $fout : $!\n";
	}

	while (<FILE>) {
		next if (/^ *\%/);
		s/outSpectrum\d\.pgfplots/$name/;
		print FOUT;
	}

	close(FOUT);
	close(FILE);
}

