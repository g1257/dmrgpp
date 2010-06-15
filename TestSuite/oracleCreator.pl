#! /usr/bin/perl -w

use strict;

runScript();

sub runScript
{
	selectOracle();
}

sub selectOracle
{
	print "Enter the number of the file to create its oracle: ";
	my $fileNum = chomp($_ =  <STDIN>);

		extractEnergy($fileNum);
#	else
#		print "The file specified cannot be found.\n";
}

sub extractEnergy
{
        my $fileAnalyzed = "data$_.txt";

        system("grep Energy $fileAnalyzed > /oracles/e$_.txt");

        print "Oracle creation is complete.\n";
}
