#!/usr/bin/perl

use strict;

=pod

PSIDOC label

PSIDOCCOPY label

=cut

my ($file) = @ARGV;
my $find = "find ../src -iname \"*.h\" -or -iname \"*.cpp\"";
defined($file) or die "USAGE: $find | $0 file\n";
my %labels;

loadLabels(\%labels);

recursiveExpand(\%labels);

replaceLabels($file,\%labels);

sub loadLabels
{
	while (<STDIN>) {
		chomp;
		my $f = $_;
		procFile(\%labels,$f);
	}
}

sub procFile
{
	my ($a,$f) = @_;
	my %labels = %$a;
	my $label = "!DISABLED";
	my $buffer = "";
	open(FILE,$f) or die "$0: Cannot open $f : $!\n";
	while (<FILE>) {
		if (/\/\* *PSIDOC *([^ ]+)/) {
			$label = $1;
			chomp($label) if ($label=~/\n$/);
			my $txt = $labels{"$label"};
			if (defined($txt)) {
				print STDERR "$0: ERROR: Label $label in file $f is duplicate\n";
				last;
			}

			next;
		}

		if (/\*\//) {
			if ($label ne "!DISABLED") {
				my @temp = ($buffer);
				$labels{"$label"} = \@temp;
			}

			$buffer = "";
			$label = "!DISABLED";
		}

		if ($label ne "!DISABLED") {
			$buffer .= $_;
		}
	}

	close(FILE);
	my $n = scalar(%labels);
	print STDERR "$0: File $f proc'ed ($n labels found so far)\n";
	%$a = %labels;
}

sub replaceLabels
{
	my ($file,$a) = @_;
	my %labels = %$a;
	my $fout = $file;
	$fout=~s/\.ptex$/\.tex/;
	die "$0: $file must have extension .ptex\n" if ($file eq $fout);

	open(FOUT,"> $fout") or die "$0: Cannot write to $fout : $!\n";
	open(FILE,$file) or die "$0: Cannot open $file : $!\n";
	while (<FILE>) {
		next if (/^[ \t]*%/);
		if (/\\ptexPaste\{([^\}]+)\}/) {
			my $label = $1;
			next if ($label eq "#1");
			my $txt = getTextFromLabel($label,$a);
			last if (!defined($txt));
			print FOUT $txt;
			next;
		}

		if (/\\ptexReadFile\{([^\}]+)\}/) {
			my $ret = open(FILE2,$1);
			if (!$ret) {
				print STDERR "$0: ERROR: Cannot read $1, line $_\n";
				last;
			}

			while (<FILE2>) {
				print FOUT;
			}

			close(FILE2);
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
		my $txt2 = expandIfNeeded($txt,$a);
		$recurse = 1 if ($txt2=~/PSIDOCCOPY/);
		my @buffer = ($txt2);
		$ptr = \@buffer;
	}

	$a = \%labels;

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
			die "Line $_\n" if (!defined($txt));
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
