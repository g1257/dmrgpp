#!/usr/bin/perl

use strict;
use warnings;
use utf8;

my ($file) = @ARGV;
defined($file) or die "USAGE: $0 filename\n";

my ($labels, $matrices) = loadLabelsAndData($file);

exit(1) if (scalar(@$labels) == 0);

my $gold_file = getGoldfilename($file);
#print STDERR "$0: Gold file is $gold_file\n";
my ($glabels, $gmatrices) = loadLabelsAndData($gold_file);
compareLabels($glabels, $labels);
compareData($gmatrices, $matrices, $glabels);

sub getGoldfilename
{
	my ($file) = @_;
	my $gold_file = `basename $file`;
	chomp($gold_file);
	my $thisname = `basename $0`;
	chomp($thisname);
	my $dir = $0;
	$dir =~ s/$thisname//;
	return "$dir/../TestSuite/gold/$gold_file";
}

sub compareData
{
	my ($gold, $steel, $labels) = @_;

	my $n = scalar(@$gold);
	my $nsteel = scalar(@$steel);
	if ($n != scalar(@$steel)) {
		die "$0: Number of gold matrices $n different than steel matrices $nsteel\n";
	}

	for (my $i = 0; $i < $n; ++$i) {
		next if (matrixEqual($gold->[$i], $steel->[$i]));
		die "$0: Matrices differ for label ".$labels->[$i]."\n";
	}
}

sub matrixEqual
{
	my ($ma, $mb) = @_;
	my $na = scalar(@$ma);
	my $nb = scalar(@$mb);
	return 0 if ($na != $nb);
	for (my $i = 0; $i < $na; ++$i) {
		next if (vectorEqual($ma->[$i], $mb->[$i]));
		return 0;
	}

	return 1;
}

sub vectorEqual
{
	my ($va, $vb) = @_;
	my $na = scalar(@$va);
	my $nb = scalar(@$vb);
	return 0 unless ($na == $nb);
	for (my $i = 0; $i < $na; ++$i) {
		next if (abs($va->[$i] - $vb->[$i]) < 1e-5);
		return 0;
	}

	return 1;
}

sub compareLabels
{
	my ($gold, $steel) = @_;
	my $n = scalar(@$gold);
	my $nsteel = scalar(@$steel);
	if ($n != scalar(@$steel)) {
		die "$0: Number of gold labels $n different than steel labels $nsteel\n";
	}

	for (my $i = 0; $i < $n; ++$i) {
		next if ($gold->[$i] eq $steel->[$i]);
		die "$0: Gold label ".$gold->[$i]." differs from ".$steel->[$i]."\n";
	}
}

sub loadLabelsAndData
{
	my ($file) = @_;
	my @labels;
	my @matrices;
	open(my $fh, "<", $file) or die "$0: Cannot open $file : $!\n";
	while (<$fh>) {
		chomp;
#PsiApp: CmdLine: /home/gonzalo/github/dmrgpp/build/src/observe -l +r -f /home/gonzalo/github/dmrgpp/src/../TestSuite/inputs/input2.ain <gs|c';c|gs>,<gs|2.0*sz;2.0*sz|gs>,<gs|n;n|gs>
		if (/CmdLine: /) {
			my @temp = split;
			my $n = scalar(@temp);
			if ($n < 6) {
				print STDERR "$0: CmdLine has error in $file\n";
				last;
			}

			my $str = $temp[$n - 1];
			$str =~ s/ //;
			@labels = split/,/, $str;
			my $nlabels = scalar(@labels);
			if ($nlabels == 0) {
				print STDERR "$0: No labels found in $file\n";
				last;
			}

			next;
		}

		my $error_flag = 0;
		foreach my $label (@labels) {
			my $line = $_;
			if ($line =~ /^$label/) {
				$_ = <$fh>;
				chomp;
				my @temp = split;
				if (scalar(@temp) != 2 or $temp[0] != $temp[1]) {
					print STDERR "$0: Malformed data for $label in $file : $. $_\n";
					$error_flag = 1;
					last;
				}

				my $ndata = $temp[0];
				my $j = 0;
				my @matrix;
				for (my $i = 0; $i < $ndata; ++$i) {
					$_ = <$fh>;
					chomp;
					my @temp = split;
					$matrix[$j++] = \@temp;
				}

				push @matrices, \@matrix;
			}
		}

		if ($error_flag) {
			@labels = ();
			last;
		}
	}

	close($fh);
	print STDERR "$0: Found ".scalar(@matrices)." correlations in $file\n";
	return (\@labels, \@matrices);
}

