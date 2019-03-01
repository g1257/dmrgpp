#!/usr/bin/perl6

use v6;

my $myself = $*PROGRAM-NAME;

sub MAIN($file)
{
	my $input = open($file, :r);
	my Int $ln = 0;

	print "##Ainur1.0\n";
	for $file.IO.lines -> $line {
		++$ln;
		my $copy = $line;
		
		$copy ~~ s/\.txt// if ($copy ~~ /OutputFile\=/);
		$copy ~~ s/\=(<-[\d\-\.]>.*$)/\=\"$0\"/;

		my $sc = ($copy ~~ /^$/) ?? "" !! ";"; 
		print "$copy$sc\n";
	}
}

