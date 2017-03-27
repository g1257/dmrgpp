#!/usr/bin/perl
#
use strict;
use warnings;

# Accepted runs between 1 minute and 1/2 hour
my ($tmin, $tmax) = (60, 1800);
# Accepted files up to 10MB
my $sizemax = 1024*1024*10; # 10 MB

my $cannotAccept = "$0: I cannot accept this run";

my ($input) = @ARGV;
defined ($input) or die "USAGE: $0 input.inp\n";

my $root = $input;

if ($root =~ /\.inp$/) {
	$root =~ s/\.inp$//;
} else {
	$input .= ".inp";
}

my $datafile = readLabel($input, "OutputFile=");
my $cout = "runFor$root.cout";
my $t = getRunTime($cout);

if ($t < $tmin || $t > $tmax) {
	print STDERR "$cannotAccept due to runtime=$t\n";
	die "$0: Runtime must be between $tmin and $tmax seconds\n";
}

my @files;

push @files, $input;

push @files, $cout;

push @files, $datafile;

my $targz = "$root.tar.gz";
tarAndGzipFiles($targz, \@files);
# Check file size

my $size = -s "$targz";
if ($size > $sizemax) {
	print STDERR "$0: $cannotAccept due to filesize=$size > maxsize=$sizemax\n";
}

sub readLabel
{
	my ($file, $label) = @_;
	my $value;
	open(FILE, "$file") or die "$0: Cannot open file $file : $!\n";
	while (<FILE>) {
		if (/^$label(.*$)/) {
			$value = $1;
			last;
		}
	}

	close(FILE);

	defined($value) or die "$0: File $file does not contain label $label\n";
	return $value;
}

sub getRunTime
{
	my ($file) = @_;
	my $v1 = readLabel($file, "UnixTimeStart=");
	my $v2 = readLabel($file, "UnixTimeEnd=");
	return ($v2 - $v1);
}

sub tarAndGzipFiles
{
	my ($targz, $files) = @_;
	if (-w "$targz") {
		my $cmd = "cp $targz $targz.bak";
		print STDERR "$0: Making backup $cmd ...\n";
		system($cmd);
	}

	my $cmd = "tar --owner=nobody --group=nobody -zcf $targz @files";
	print STDERR "$0: Running $cmd ...\n";
	my $ret = system($cmd);
	if ($ret != 0) {
		die "$0: tar command failed\n";
	}

	print STDERR "$0: File $targz created successfully\n";
}
 
	

	

