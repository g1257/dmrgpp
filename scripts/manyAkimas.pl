#!/usr/bin/perl -w

use strict;

#dmrgpp52.o165
my @files=@ARGV;
buildFilenames() if (!defined($files[0]));
print STDERR "Found $#files\n";
my $firstSite = 0;
my $n = 6;

foreach my $file (@files) {
	for (my $site=$firstSite;$site<$n;$site++) {
		my $fout = $file;
		$fout =~ s/\..*$//;
		$fout = $fout.".nd".$site;
		die "$0: input and output files are the same\n" if ($fout eq $file);
		system("perl getTimeObs.pl $site $file  > $fout");
	}
}

print STDERR "Sites done successfully\n";

foreach my $file (@files) {
	my $fin = $file;
	$fin =~ s/\..*$/\.nd/;
	my $fout = $fin."t";
	die "$0: input and output files are the same\n" if ($fout eq $file);
	my $x = $n;
	print STDERR  "$0: Procing $file with x = $x...\n";
	system("perl akimaSumBatch.pl $fin $x \"\" > $fout");
}

sub buildFilenames
{
	open(PIPE,"ls dmrg14*.o* |") or die "$0: Cannot open pipe: $!\n";
	my $counter = 0;
	while(<PIPE>) {
		chomp;
		$files[$counter++] = $_ if (-r $_);
	}
	close(PIPE);

}

