#!/usr/bin/perl -w

use strict;
use warnings;
use File::Find::Rule;

my ($directory)=@ARGV;

defined($directory) or die "$0: USAGE: $0 directory\n";

my @files = File::Find::Rule->file()->name("*.cpp", "*.h")
                                ->in($directory);
my $counter=0;
my $fileCounter = 0;
foreach my $file (@files) {
	$counter += countThisFile($file);
	$fileCounter++;
}

print "Total Lines=$counter\n";
print "Total Files=$fileCounter\n";

sub countThisFile
{
	my ($file)=@_;
	my $s=0;
	#my $printEach = 1000;
	my $multiLineComments=0;
	open(FILE,$file) or die "Cannot open file $file : $!\n";
	while(<FILE>) {
		chomp;
		s/\t//g;
		s/ //g;
		next if ($_ eq "");
		next if (/^\/\//);
		next if (/^\/\*.*\*\/$/);
		$multiLineComments=1 if (/^\/\*/);
		$multiLineComments=0 if (/\*\/$/);
		next if ($multiLineComments);
		$s++;
		#print "*${file}*::$_\n" if (($s % $printEach)==0);
	}
	close(FILE);
	return $s;
}



