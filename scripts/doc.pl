#!/usr/bin/perl

use strict;

=pod

PSIDOC label

PSIDOCCOPY label

=cut

my ($file, $doc_root, $out_dir) = @ARGV;
my $find = "find ../src -iname \"*.h\" -or -iname \"*.cpp\"";
defined($file) or die "USAGE: $find | $0 file\n";
my %labels;

loadFiles(\%labels, $file);

my @lines;
loadLines(\@lines);

loadLabels(\%labels, \@lines);

recursiveExpand(\%labels);

replaceLabels($file, \%labels);

sub loadLines
{
	my ($lines) = @_;
	while (<STDIN>) {
		chomp;
		my $f = $_;
		loadLinesThisFile($lines, $f);
	}
}

sub loadFiles
{
	my ($a, $f) = @_;
	my %labels = %$a;

	open(FILE, "<", $f) or die "$0: Cannot open $f : $!\n";

	while (<FILE>) {
	    if (/\\ptexReadFile\{([^\}]+)\}/) {
		        my $file = $1;
			if ($doc_root != "") {
			    $file = $doc_root . "/" . $file;
			}
			my $ret = open(FILE2, "<", $file);
			if (!$ret) {
				close(FILE);
				die "$0: ERROR: Cannot read $file, line $_\n";
			}

			my $isMdFile = ($file=~/\.md$/);
			my $buffer = "";
			while (<FILE2>) {
				if ($isMdFile) {
					s/^# +(.*)/\\chapter\{$1\}/;
					s/^## +(.*)/\\section\{$1\}/;
					s/^### +(.*)\n$/\\subsection\{$1\}/;
					s/\\\@/\@/g;
					s/<b>/\{\\bf /;
					s/<\/b>/\}/;

					s/<pre>/\\begin\{verbatim\}\n/;
					s/<\/pre>/\\end\{verbatim\}\n/;
					s/<code>/\{\\tt /g;
					s/<\/code>/\}/g;
				}

				$buffer .= $_;
			}

			close(FILE2);

			my $label = getLabelForFile($file);
			my @temp = ($buffer);
			$labels{"$label"} = \@temp;
			next;
		}
	}

	close(FILE);
	%$a = %labels;
}

sub loadLinesThisFile
{
	my ($lines, $f) = @_;

	open(FILE, "<", $f) or die "$0: Cannot open $f : $!\n";
	while (<FILE>) {
		chomp;
		push @$lines, $_;
	}

	close(FILE);
}

sub loadLabels
{
	my ($a, $lines) = @_;
	my %labels = %$a;
	my $label = "!DISABLED";
	my $additional = "";
	my $buffer = "";
	my $nlines = scalar(@$lines);

	for (my $i = 0; $i < $nlines; ++$i) {
		$_ = $lines->[$i];
		if (/\/\* *PSIDOC +(.+$)/) {
			my $rest = $1;
			chomp($rest);
			($label, $additional) = procPsidocName($rest);
			my $txt = $labels{"$label"};
			if (defined($txt)) {
				die "$0: ERROR: Label $label is duplicate\n";
			}

			next;
		}

		if (/\*\//) {
			if ($label ne "!DISABLED") {
				my $inlabel = $label."::";
				$buffer =~ s/PSIDOCCOPY \$/PSIDOCCOPY ${inlabel}/g;
				my @temp = ($buffer);
				$labels{"$label"} = \@temp;
				my $proto = captureFirstProtoBelow($i + 1, \@lines);
				if ($proto) {
					my @temp = ($proto);
					my $name = $inlabel."FirstProtoBelow";
					$labels{"$name"} = \@temp;
				}
			}

			$buffer = "";
			$label = "!DISABLED";
			$additional = "";
		}

		if ($label ne "!DISABLED") {
			$buffer .= $_."\n";
		}
	}

	my $n = scalar(%labels);
	print STDERR "$0: $n labels found\n";

	%$a = %labels;
}

sub procPsidocName
{
	my ($nameLike) = @_;
	my @temp = split/[ \t]/,$nameLike;
	my $n = scalar(@temp);
	die "$0: procPsidocName empty\n" if ($n == 0);
	return $nameLike if ($n == 1);
	die "$0: procPsidocName more than one additional\n" if ($n > 2);
	return @temp;
}

sub captureFirstProtoBelow
{
	my ($ind, $lines) = @_;
	my $nlines = scalar(@$lines);
	my $buffer = "";
	for (my $i = $ind; $i < $nlines; ++$i) {
		$_ = $lines->[$i];
		last if (/\/\*/);
		next if (/^ *\/\//);
		$buffer .= $_."\n";
		last if (/\;/);
	}

	return $buffer if ($buffer eq "");

	$buffer =~ s/\t/  /g;
	$buffer = "\\begin{lstlisting}\n$buffer\n";
	$buffer .= "\\end{lstlisting}\n";

	return $buffer;
}

sub replaceLabels
{
	my ($file, $a) = @_;
	my $fout = $file;
	if ($out_dir != "") {
	    $fout =~ s/.*\///; #chop the beginning of the path
	    $fout =  $out_dir . "/" . $fout;
	}
	$fout=~s/\.ptex$/\.tex/;
	die "$0: $file must have extension .ptex\n" if ($file eq $fout);

	open(FOUT, ">", "$fout") or die "$0: Cannot write to $fout : $!\n";
	open(FILE, "<", $file) or die "$0: Cannot open $file : $!\n";
	while (<FILE>) {
		next if (/^[ \t]*%/);
		if (/\\ptexPaste\{([^\}]+)\}/) {
			my $label = $1;
			next if ($label eq "#1");
			my $txt = getTextFromLabel($label,$a);
			if (!defined($txt)) {
				$txt = labelNotFound($label);
			}

			print FOUT $txt;
			next;
		}

		if (/\\ptexReadFile\{([^\}]+)\}/) {
			my $file=$1;
			my $label = getLabelForFile($file);
			my $txt = getTextFromLabel($label,$a);
			if (!defined($txt)) {
				$txt = labelNotFound($label,$file);
			}

			print FOUT $txt;
			next;
		}

		print FOUT;
	}


	close(FILE);
	close(FOUT);

	print STDERR "$0: File $fout written\n";
}

sub recursiveExpand
{
	my ($a) = @_;
	my %labels = %$a;

	my $recurse = 0;
	foreach my $key (keys %labels) {
		my $ptr = $labels{"$key"};
		defined($ptr) or die "$0: Label $key has error\n";
		scalar(@$ptr) > 0 or die "$0: Label $key has error\n";
		my $txt = $ptr->[0];
		next unless ($txt=~/PSIDOCCOPY/);
		my $txt2 = expandIfNeeded($txt,$a);
		$recurse = 1 if ($txt2=~/PSIDOCCOPY/);
		my @buffer = ($txt2);
		$ptr = \@buffer;
		$labels{"$key"}=$ptr;
	}

	%$a = %labels;

	recursiveExpand($a) if ($recurse);
}

sub expandIfNeeded
{
	my ($txt,$a) = @_;
	my %labels = %$a;
	my @temp = split/\n/,$txt;

	my $n = scalar(@temp);
	my $buffer = "";
	for (my $i = 0; $i < $n; ++$i) {
		$_ = $temp[$i];
		if (/PSIDOCCOPY +([^ ]+)/) {
			my $label = $1;
			chomp($label) if ($label=~/\n$/);
			my $txt = getTextFromLabel($label,$a);
			if (!defined($txt)) {
				print STDERR "$0: ERROR: line $_, $label not found\n";
				$txt = notFoundLabel($label);
			}

			$buffer .= $txt;
			next;
		}

		$buffer .= $_."\n";
	}

	return $buffer;
}

sub getTextFromLabel
{
	my ($label,$a)=@_;
	my %labels = %$a;
	my $ptr = $labels{"$label"};
	my $txt;
	if (!defined($ptr) || scalar(@$ptr) < 1) {
		print STDERR "$0: ERROR: Label $label not found\n";
		return $txt;
	}

	return $ptr->[0];
}

sub getLabelForFile
{
	my ($file) = @_;

	$file=~s/\./DOT/g;
	$file=~s/\//SLASH/g;
	return "FILE$file";
}

sub labelNotFound
{
	my ($label,$file) = @_;
	my $str = (defined($file)) ? " in $file " : "";
	return "{\\bf\\textcolor{red}{ERROR: Label not found$str}}\n";
}

