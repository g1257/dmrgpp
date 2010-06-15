#! /usr/bin/perl -w

use strict;
use Getopt::Long;
my ($inputFile,$inputDir,$outputFile,$extension) = (undef,undef,"hashTable.txt","spec");

GetOptions("f=s" => \$inputFile, "d=s" => \$inputDir, "o=s" => \$outputFile, "ext=s" => \$extension);

if(!defined($inputFile) && !defined($inputDir)) {
	print "\nERROR: An input file or directory is needed!\n";
} elsif(defined($inputFile) && defined($inputDir)) {
	print "\nERROR: Only one type of input can be specified!\n";
} else {
	
	if (defined($inputDir)) {
		hashDirectory($inputDir,$outputFile) if(validateDirectory($inputDir));
	} else {
		hashFile($inputFile,$outputFile) if(validateFile($inputFile));
	}
}

sub validateDirectory
{
	my ($dir) = @_;
	
	if($dir !~ /\/$/) {
		$inputDir = "$dir/";
	}
	
	if(-d $inputDir) {
		return 1;
	} else {
		print "ERROR: The directory does not exists.\n";
		return 0;
	}
}

sub validateFile
{
	my ($file) = @_;

	if(-e $file) {
		return 1;
	} else {
		print "ERROR: The input file does not exists.\n";
		return 0;
	}
}

sub hashDirectory
{
	my ($dir,$output) = @_;
	my $iterator;
	my @searchFiles = glob "$dir*.$extension";
	
	foreach $iterator(@searchFiles) {
		system("md5sum $iterator >> $output");
	}
	
	print "Hashes were created...\n";
}

sub hashFile
{
	my ($file,$output) = @_;
	
	my $h = system("md5sum $file"); # >> $output");
#	my @temp = split(/ /,$_);
#	my $hashKey = $temp[0];
	open($output,">>$h");
	print "Hash was created...\n";
}

