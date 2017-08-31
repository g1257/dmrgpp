#!/usr/bin/perl -w

use strict;

my ($file)=@ARGV;
my $n = 12;
my $offset = 0.2;
my @times;
my %val;
my @nd;

my $sd = getSuperDensity($file);

#Note: the two extreme are not available since
# DMRG++ cannot compute them
my $firstSite = 1;
for (my $site=$firstSite;$site<$n-1;$site++)
{
	doOneSite($site);
}

printTotalNd($file,$sd);
printAllSites($file,$sd);

sub getOutputName
{
	my ($file,$ext) = @_;
	my $fout = $file;
	$fout =~ s/\..*$/\.${ext}/;
	die "$0: Refusing to overwrite $fout\n" if ($fout eq $file);
	#die "$0: Refusing to overwrite $fout\n" if (-e $fout);
	return $fout;
}

sub printTotalNd
{
	my ($file,$sd)=@_;
	my $fout = getOutputName($file,"ndt");

	open(FOUT, ">", "$fout") or die
		"Cannot open $fout for writing: $!\n";

	foreach my $t (@times) {
		$_=$val{$t}/$sd;
		print FOUT "$t $_\n";
	}
	close(FOUT);
}

sub printAllSites
{
	my ($file,$sd)=@_;
	my $fout = getOutputName($file,"nds");
	open(FOUT, ">", "$fout") or die
		"Cannot open $fout for writing: $!\n";

	# first row (labels:):
	print FOUT "#times/sites ";
	foreach my $t (@times) {
		print FOUT "$t ";
	}
	print FOUT "\n";
	
	foreach my $t (@times) {
		print FOUT "$t ";
		for (my $site=$firstSite;$site<$n-1;$site++) {
			$_=$nd[$site]{$t}/$sd; # + $offset * $site;
			print FOUT "$_ ";
		}
		print FOUT "\n";
	}
	close(FOUT);

}

sub getSuperDensity
{
	my ($site)=@_;
	my $sd;
	open(FILE, "<", $file) or die "Cannot open file $file: $!\n";
	while(<FILE>) {
		if (/SuperDensity.*=\(([^,]+),/) {
			$sd = $1;
			last;
		}
	}
	close(FILE);
	defined $sd or die "SuperDensity is not defined\n";
	return $sd;
}

sub doOneSite
{
	my ($site)=@_;

	system("perl getTimeObs.pl $site < $file > out1");
	system("perl gatherTimes.pl < out1 > out ");
	open(FILE, "<", "out") or die "Cannot open file out: $!\n";
	my $counter = 0;
	my @timesNow;
	while(<FILE>) {
		my @temp=split;
		$timesNow[$counter++]=$temp[0];
		if (defined($val{$temp[0]})) {
			$val{$temp[0]} += $temp[1];
		} else {
			$val{$temp[0]} = $temp[1];
		}
		$nd[$site]{$temp[0]}=$temp[1];
	}
	close(FILE);
	if ($site==$firstSite) {
		@times = @timesNow;
		print STDERR "Setting $#times\n";
	} else {
		updateTimes(\@timesNow,$#timesNow+1,$site);
	}
}

sub updateTimes
{
	my ($timesNow,$n,$site)=@_;
	# If there's a time in @times that's not in @timesNow
	# then delete it
	my $counter=0;
	my @newTimes;
	foreach my $t (@times) {
		if (isInVector($timesNow,$n,$t)) {
			$newTimes[$counter++]=$t;
		} else {
			#print STDERR "Eliminated $t when site=$site\n";
			
			#if ($t==0.6) {
			#	die "timesnow: @$timesNow\n";
			#}
		}
	}
	@times = @newTimes;
	print STDERR "Updated to $#times\n";
}

sub isInVector
{
	my ($v,$n,$what)=@_;
	for (my $i=0;$i<$n;$i++) {
		return 1 if ($v->[$i]==$what);
	}
	return 0;
}


