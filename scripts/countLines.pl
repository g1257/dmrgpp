#!/usr/bin/perl -w

use strict;

open(PIPE,"find ../src -iname \"*.h\" -o -iname \"*.cpp\" | ") or die "Cannot open pipe: $!\n";

my $counter=0;
my $fileCounter = 0;
while(<PIPE>) {
	next if ($_ eq "" or $_ eq "\n");
	$counter += countThisFile($_);
	$fileCounter++;
}

close(PIPE);
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



