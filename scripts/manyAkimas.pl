#!/usr/bin/perl -w

use strict;

#dmrgpp52.o165
my @files=@ARGV;
buildFilenames() if (!defined($files[0]));

my $firstSite = 0;
my $n = 6;
my $whatState = "time";
my $whatObservable = "nd";

foreach my $file (@files) {
	for (my $site=$firstSite;$site<$n;$site++) {
		my $fout = $file;
		$fout =~ s/\..*$//;
		$fout = $fout.".${whatState}-${whatObservable}".$site;
		die "$0: input and output files are the same\n" if ($fout eq $file);
		system("perl getTimeObs.pl $site $file $whatState $whatObservable> $fout");
	}
}

foreach my $file (@files) {
	my $fin = $file;
	$fin =~ s/\..*$/\.${whatState}-${whatObservable}/;
	my $fout = $fin."t";
	die "$0: input and output files are the same\n" if ($fout eq $file);
	my $x = $n;
	system("perl akimaSumBatch.pl $fin $x \"\" > $fout");
}

sub buildFilenames
{
	open(PIPE,"ls dmrgpp*.o* |") or die "Cannot open pipe: $!\n";
	$_=<PIPE>;
	@files = split;
	close(PIPE);

}

