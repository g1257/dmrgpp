#!/usr/bin/perl
#
use strict;
use warnings;
use utf8;

my ($uroot, $runLatex, $lOrUp) = @ARGV;

defined($uroot) or die "USAGE: $0 root [runLatex=1] [L | U | sz]\n";
defined($runLatex) or $runLatex = 1;

if (!defined($lOrUp)) {
	if ($uroot =~ /(^[^\.]+)\.pgfplots/) {
		$uroot = $1;
	}

	my $file = "sample.tex";
	my $fout = "$uroot.tex";
	my $name = "$uroot.pgfplots";
	fromTexToTex($fout, $file, $name, $runLatex);
	exit(0);
}

my $root = "$uroot$lOrUp"."_ky";

doFile(0, $root, $runLatex);
doFile(1, $root, $runLatex);

sub doFile
{
	my ($ind, $root, $runLatex) = @_;

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
	fromTexToTex($fout, $file, $name, $runLatex);
}

sub fromTexToTex
{
	my ($fout, $file, $name, $runLatex) = @_;

	my $dirForTex = $0;
	$dirForTex =~ s/pgfplot\.pl$//;
	die "$0: Not a directory $dirForTex\n" unless (-d "$dirForTex");

	system("cp $dirForTex/sample.tex .") unless (-r "sample.tex");
	system("cp $dirForTex/palette.tex .") unless (-r "palette.tex");

	copyAndEdit($fout, $file, $name);

	return if (!$runLatex);

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

