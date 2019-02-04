#!/usr/bin/perl
#
use strict;
use warnings;
use utf8;

my ($dir, $h, $u, $lOrUp) = @ARGV;
defined($lOrUp) or die "USAGE: $0 dir h u L | U\n";

die "$0: $dir not a dir\n" unless (-d "$dir");

my $root = "16x4u$u"."nih$h$lOrUp"."akwladder0ky";

doFile(0, $root, $dir);
doFile(1, $root, $dir);

sub doFile
{
	my ($ind, $root, $dir) = @_;

	my $file = "$dir/outSpectrum$ind.pgfplots";
	my $fout = "$root$ind.pgfplots";
	unlink("$fout");
	my $cmd = "cp $file $fout";
	die "$0: Failed to create $fout\n" if (-r "$fout");
	system("$cmd");

	$file = "sample.tex";
	$fout = "$root$ind.tex";
	copyAndEdit($fout, $file, $ind);

	$cmd = "pdflatex $fout";
	system("$cmd");
}

sub copyAndEdit
{
	my ($fout, $file, $ind) = @_;
	open(FILE, "<", "$file") or die "$0: Cannot open $file : $!\n";
	my $ret = open(FOUT, ">", "$fout");
	if (!$ret) {
		close(FILE);
		die "$0: Cannot open write to $fout : $!\n";
	}

	my $name = $fout;
	$name =~ s/tex$/pgfplots/;
	while (<FILE>) {
		next if (/^ *\%/);
		s/outSpectrum\d\.pgfplots/$name/;
		print FOUT;
	}

	close(FOUT);
	close(FILE);
}

