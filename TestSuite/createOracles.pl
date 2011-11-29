#!/usr/bin/perl -w

my @roots=("e","timeEvolution","operatorC","operatorN","operatorSz");

my ($doIt)=@ARGV;
$doIt = 0 if (!defined($doIt));

foreach (@roots) {
	doThisRoot($_);
}

sub doThisRoot
{
	my ($root)=@_;
	for ($i=0;$i<800;$i++) {
		my $f = "$root$i.txt";
		next unless (-r "results/$f");
		system("cp results/$f oracles/$f") if ($doIt);
		print STDERR "cp results/$f oracles/$f\n";
	}
}

