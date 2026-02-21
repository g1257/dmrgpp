#!/usr/bin/perl

use warnings;
use strict;
use utf8;

my ($cmakelists) = @ARGV;
defined($cmakelists) or die "USAGE: $0 cmakelists.txt\n";
(-r "$cmakelists") or die "$0: File $cmakelists not readable\n";

my @files = buildFileList();

my $h = extractEnergies(\@files);

writeCmakeLists($cmakelists, $h);

# add_test(NAME input4 COMMAND ./dmrg -f ${PATH_TO_INPUTS}/input4.ain)
sub writeCmakeLists
{
	my ($cmakelists, $h) = @_;
	open(my $fh, "<", $cmakelists) or die "$0: Cannot open $cmakelists : $!\n";
	while (<$fh>) {
		chomp;
		my $line = $_;
		if (/\.\/dmrg / and /(^ *add_test *\(.+)(input\d+\.ain) *([^\)]*)\) *$/) {
			my $pre = $1;
			my $file = $2;
			my $post = $3;
			my $cout = buildCout($file);
			my $hash = $h->{$cout};
			my ($what, $value) = ($hash->{type}, $hash->{value});
			if ($what) {
				$line = "$pre$file";
				$line .= " $post" if ($post);
				$line .= ") # $what";
				if ($value) {
					$line .= " $value";
				}
			}
		}

		print "$line\n";
	}

	close($fh);
}

sub buildCout
{
	my ($input) = @_;
	$input =~ s/\.ain$//;
	return "runFor$input.cout";
}

sub buildFileList
{
	my @files;
	my $cmd = "find src -iname \"runForinput*.cout\"";
	open(PIPE, "$cmd |") or die "$0: Cannot open pipe $cmd : $!\n";
	while (<PIPE>) {
		chomp;
		my $file = $_;
		next unless $file =~ /runForinput\d+\.cout/;
		push @files, $file;
	}

	close(PIPE);
	print STDERR "$0: Found ".scalar(@files)." files\n";
	return @files;
}

sub extractEnergies
{
	my ($files) = @_;
	my $h;
	foreach my $file (@files) {
		my ($e, $d) = readLabel($file, "Found lowest eigenvalue= ");
		my $name = `basename $file`;
		chomp($name);
		if (!defined($e)) {
			$h->{$name} = {type => "none"};
		} else {
			$h->{$name} = {type => "energy", value => $e};
		}
	}

	return $h;
}

sub readLabel
{
	my ($file, $label) = @_;
	open(my $fh, "<", $file) or die "$0: Cannot open $file : $!\n";
	my ($e, $d);
	my $eprev;
	while (<$fh>) {
		chomp;
		if (/$label([^ ]+)/) {
			$e = $1;
			$d = abs($eprev - $e) if (defined($eprev));
			$eprev = $e;
		}
	}

	close($fh);
	return ($e, $d);
}
