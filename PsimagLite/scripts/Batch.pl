#!/usr/bin/perl
#PBS -N mpi
#PBS -q batch
#PBS -l nodes=2:ppn=16
#PBS -l walltime=1:00:00
#PBS -l vmem=90gb

use strict;
use warnings;
use utf8;

chdir $ENV{"PBS_O_WORKDIR"};

#date
my $np=4;
my $nodefile = $ENV{"PBS_NODEFILE"};
my $h = getHostString($nodefile, $np);

#my $realCmd = "hostname";
my $realCmd = "./internode 4";

my $command = "LD_LIBRARY_PATH=\${LD_LIBRARY_PATH}:/usr/lib64/openmpi/lib ";
$command .= "/usr/lib64/openmpi/bin/mpiexec -n $np -host \"$h\"  $realCmd";
system("$command");

sub getHostString
{
	my ($file, $mpiJobs) = @_;

	my %h;
	my $firstNode;
	open(FILE, "<", "$file") or die "$0: Cannot open $file : $!\n";
	while (<FILE>) {
		chomp;
		next if (/^#/ or $_ eq "");
		my $f = $_;
		$firstNode = $f;
		if (!defined($h{"$f"})) {
			$h{"$f"} = 1;
		} else {
			++$h{"$f"};
		}
	}

	close(FILE);

	my $nodes = scalar(keys %h);
	die "$0: Nodes $nodes must be a multiple of mpiJobs $mpiJobs\n" if ($mpiJobs % $nodes != 0);

	die "$0: No nodes!\n" if ($nodes == 0 or !defined($firstNode));
	my $ppn = $h{"$firstNode"};
	my $repeat = $mpiJobs/$nodes;
	return getHostString2($repeat, \%h);
}

sub getHostString2
{
	my ($repeat, $h) = @_;
	my $firstcall = 1;
	my $str = "";
	for (my $i = 0; $i < $repeat; ++$i) {
		foreach my $key (sort keys %$h) {
			my $n = $h->{"$key"};
			next if ($n == 0);
			if ($firstcall) {
				$firstcall = 0;
			} else {
				$str .= ",";
			}

			$str .= "$key";
		}
	}

	return $str;
}

