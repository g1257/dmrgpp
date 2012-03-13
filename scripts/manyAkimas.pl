#!/usr/bin/perl -w

use strict;
use lib '/home/gz1/programs/github/dmrgpp/scripts';
use GetTimeObs;
use AkimaSumBatch;

my $whatState = shift;
my $whatObservable = shift;
my $start = shift;
my $end = shift;
my $points = shift;
my $n;

my @files = @ARGV;
buildFilenames() if (!defined($files[0]));

foreach my $file (@files) {
	my @allSites = getAllSites($file);
	@allSites = sort {$a <=> $b} @allSites;
	my $firstSite = $allSites[0];
	$n = $allSites[$#allSites]  + 1;
	print STDERR "Lowest site =  $firstSite , largest site = $n\n";
	for (my $site=$firstSite;$site<$n;$site++) {
		my $fout = $file;
		$fout =~ s/\..*$//;
		$fout = $fout.".${whatState}-${whatObservable}".$site;
		die "$0: input and output files are the same\n" if ($fout eq $file);
		my $fh;
		open($fh,"> $fout") or die "Cannot open file $fout for writing: $!\n";
		GetTimeObs::main($fh,$site,$file,$whatState,$whatObservable);
		close($fh);
	}
}

foreach my $file (@files) {
	my $fin = $file;
	$fin =~ s/\..*$/\.${whatState}-${whatObservable}/;
	my $fout = $fin."t";
	die "$0: input and output files are the same\n" if ($fout eq $file);
	my $x = $n;
	my $fh;
	open($fh,"> $fout") or die "Cannot open file $fout for writing: $!\n";
	AkimaSumBatch::main($fh,$start,$end,$points,$fin,$x,"");
	close($fh);
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

