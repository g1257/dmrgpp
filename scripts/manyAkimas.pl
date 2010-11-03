#!/usr/bin/perl -w

use strict;

#dmrgpp52.o165
my @files;
buildFilenames();

my $firstSite = 1;
my $n = 12;

foreach my $file (@files) {
	for (my $site=$firstSite;$site<$n-1;$site++) {
		my $fout = $file;
		$fout =~ s/\..*$//;
		$fout = $fout.".nd".$site;
		die "$0: input and output files are the same\n" if ($fout eq $file);
		system("perl getTimeObs.pl $site $file  > $fout");
	}
}

foreach my $file (@files) {
	my $fin = $file;
	$fin =~ s/\..*$/\.nd/;
	my $fout = $fin."t";
	die "$0: input and output files are the same\n" if ($fout eq $file);
	my $x = $n-2;
	system("perl akimaSumBatch.pl $fin $x \"\" > $fout");
}

sub buildFilenames
{
	open(PIPE,"ls dmrgpp*.o* |") or die "Cannot open pipe: $!\n";
	$_=<PIPE>;
	@files = split;
	close(PIPE);

}

