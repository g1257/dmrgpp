#!/usr/bin/perl -w

use strict;

#dmrgpp52.o165

my $whatState = shift;
my $whatObservable = shift;
my $n;

my @files = @ARGV;
buildFilenames() if (!defined($files[0]));

foreach my $file (@files) {
	my @allSites = getAllSites($file);
	@allSites = sort @allSites;
	my $firstSite = $allSites[0];
	$n = $allSites[$#allSites] - $firstSite + 1;
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

sub getAllSites
{
	my ($file) = @_;
	open(FILE,$file) or die "Cannot open $file: $!\n";
	my $counter = 0;
	my @a;
	while(<FILE>) {
		next unless (/^\d+/);
		my @temp = split;
		$a[$counter++]=$temp[0];
	}
	close(FILE);
	#print STDERR "@a\n";
	#die "Testing\n";
	return @a;
}

