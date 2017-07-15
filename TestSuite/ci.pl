#!/usr/bin/perl

use strict;
use warnings;
use Getopt::Long qw(:config no_ignore_case);
use Ci;

my ($valgrind,$workdir,
    $restart,$ranges,$postprocess,$noSu2,$help);
my %submit;
GetOptions(
'S=s' => \$submit{"command"},
'delay=s' => \$submit{"delay"},
'valgrind=s' => \$valgrind,
'w=s' => \$workdir,
'R' => \$restart,
'P' => \$postprocess,
'n=s' => \$ranges,
'nosu2' => \$noSu2,
'h' => \$help) or die "$0: Error in command line args, run with -h to display help\n";

if (defined($help)) {
	print "USAGE: $0 [options]\n";
	print "\tIf no option is given creates inputs and batches for all ";
	print "no-restart tests\n";
	print "\t-S command\n";
	print "\t\tAfter creating inputs and batches, submit them to the queue\n";
	print "\t\tusing batchDollarized.pbs as template, and\n";
	print "\t\twith command command, usually command is qsub but you can also use\n";
	print "\t\tbash to run in the command line without a batching system.\n";
	print "\t-delay delay\n";
	print "\t\tDelay in seconds between subsequent submissions.\n";
	print "\t-n n\n";
	print Ci::helpFor("-n");
	print Ci::helpFor("-w");
	print "\t-R\n";
	print "\t\tRun restart tests only\n";
	print "\t-P\n";
	print "\t\tRun postprocess only\n";
	print "\t--valgrind tool\n";
	print "\t\tRun with valgrind using tool tool\n";
	print Ci::helpFor("-nosu2");
	print Ci::helpFor("-h");
	exit(0);
}

defined($submit{"command"}) or $submit{"command"} = "";
defined($noSu2) or $noSu2 = 0;

defined($valgrind) or $valgrind = "";
defined($workdir) or $workdir = "tests";
defined($restart) or $restart = 0;
defined($postprocess) or $postprocess = 0;

my $templateBatch = "batchDollarized.pbs";
my @tests = Ci::getTests();

prepareDir();

my $total = scalar(@tests);

my @inRange = Ci::procRanges($ranges, $total);
my $rangesTotal = scalar(@inRange);

die "$0: No tests specified under -n\n" if ($rangesTotal == 0);

for (my $j = 0; $j < $rangesTotal; ++$j) {
	my $i = $inRange[$j];
	die "$0: out of range $i >= $total\n" if ($i >= $total);
	my $n = $tests[$i];

	my $ir = isRestart("../inputs/input$n.inp",$n);
	if ($restart and !$ir) {
		print STDERR "$0: WARNING: Ignored test $n ";
		print STDERR "because it's a restart test and ";
		print STDERR "you did NOT specify -R\n";
		next;
	}

	if (!$restart and $ir) {
		print STDERR "$0: WARNING: Ignored test $n ";
		print STDERR "because it's NOT a restart test and ";
		print STDERR "you specified -R\n";
		next;
	}

	my $isSu2 = Ci::isSu2("../inputs/input$n.inp",$n);
	if ($isSu2 and $noSu2) {
		print STDERR "$0: WARNING: Ignored test $n ";
		print STDERR "because it's NOT an SU(2) test and ";
		print STDERR "you specified -nosu2\n";
		next;
	}

	my %ciAnnotations = getCiAnnotations("../inputs/input$n.inp",$n);
	my $whatObserve = $ciAnnotations{"observe"};
	my $whatObserveN = scalar(@$whatObserve);
	if ($whatObserveN > 0) {
		print "|$n| has $whatObserveN observe lines\n";
		if ($postprocess) {
			runObserve($n,$whatObserve,\%submit);
			next;
		}
	}

	my $whatDmrg = $ciAnnotations{"dmrg"};
	my $extraCmdArgs = findArguments($whatDmrg);
	procTest($n,$valgrind,\%submit,$extraCmdArgs);
}

sub findArguments
{
	my ($a) = @_;
	my $n = scalar(@$a);
	for (my $i = 0; $i < $n; ++$i) {
		if ($a->[$i] =~/^arguments=(.+$)/) {
			return $1;
		}

		die "$0: #ci dmrg annotation: $a->[$i] not understood\n";
	}

	return "";
}

sub runObserve
{
	my ($n,$what,$submit) = @_;
	my $whatN = scalar(@$what);
	my $cmd = "";
	for (my $i = 0; $i < $whatN; ++$i) {
		$cmd .= runObserveOne($n,$i,$what->[$i],$submit);
		$cmd .= "\n";
	}

	my $batch = createBatch($n,$cmd);
	submitBatch($submit, $batch) if ($submit->{"command"} ne "");
}

sub runObserveOne
{
	my ($n,$ind,$what,$submit) = @_;
	# what == arguments=something
	my $args;
	if ($what =~ /^arguments=(.+$)/) {
		$args = $1;
	}

	defined($args) or die "$0: observe must have arguments\n";

	my $output = "observe$n.txt";
	unlink($output) if ($ind == 0);
	my $cmd = "./observe -f ../inputs/input$n.inp \"$args\" >> $output";
	print "|$n|: postTest $cmd\n";
	return $cmd;
}

sub getCiAnnotations
{
	my ($file,$n) = @_;
	open(FILE, "$file") or return "";
	my @whatObserve;
	my @whatDmrg;
	my $counter = 0;
	while (<FILE>) {
		chomp;
		if (/^\#ci observe (.*$)/) {
			push (@whatObserve, "$1");
			next;
		}

		if (/^#ci dmrg (.*$)/) {
			push (@whatDmrg, "$1");
			next;
		}
	}

	close(FILE);
	my %h;
	$h{"dmrg"} = \@whatDmrg;
	$h{"observe"} = \@whatObserve;
	return %h;
}

sub isRestart
{
	my ($file,$n) = @_;
	open(FILE, "$file") or return 0;
	my $so;
	while (<FILE>) {
		chomp;
		if (/SolverOptions=(.*$)/) {
			$so = $1;
			last;
		}
	}

	close(FILE);
	defined($so) or return 0;
	if ($so =~/restart/) { return 1; }
	return 0;
}

sub procTest
{
	my ($n,$tool,$submit,$extraCmdArgs) = @_;
	my $valgrind = ($tool eq "") ? "" : "valgrind --tool=$tool ";
	$valgrind .= " --callgrind-out-file=callgrind$n.out " if ($tool eq "callgrind");
	my $cmd = "$valgrind./dmrg -f ../inputs/input$n.inp $extraCmdArgs &> output$n.txt";
	my $batch = createBatch($n,$cmd);
	submitBatch($submit, $batch) if ($submit->{"command"} ne "");
}

sub createBatch
{
	my ($ind,$cmd) = @_;
	my $file = "Batch$ind.pbs";
	open(FOUT,">$file") or die "$0: Cannot write to $file: $!\n";

	open(FILE,"../$templateBatch") or die "$0: Cannot open ../$templateBatch: $!\n";

	while (<FILE>) {
		while (/\$\$([a-zA-Z0-9\[\]]+)/) {
			my $line = $_;
			my $name = $1;
			my $str = "\$".$name;
			my $val = eval "$str";
			defined($val) or die "$0: Undefined substitution for $name\n";
			$line =~ s/\$\$$name/$val/;
			$_ = $line;
		}

		s/cd +\$PBS_O_WORKDIR//;
		print FOUT;
	}

	close(FILE);
	close(FOUT);

	print STDERR "$0: $file written\n";
	return $file;
}

sub submitBatch
{
	my ($submit, $batch) = @_;
	my $extra = $submit->{"extra"};
	my $delay = $submit->{"delay"};
	defined($delay) or die "$0: Please say -delay delay, where delay > 0 in seconds\n";
	($delay > 0) or die "$0: delay must be positive\n";
	defined($extra) or $extra = "";
	sleep(1);
	print STDERR "$0: Submitted $batch $extra $batch\n";

	my $qsub = $submit->{"command"};
	($qsub eq "bash" or $qsub eq "qsub") or die "$0: -S qsub | bash\n";
	my $ret = `$qsub $extra $batch`;
	defined($ret) or $ret = "";
	chomp($ret);
	sleep($delay);
	return $ret;
}

sub prepareDir
{
	my $b = (-r "$workdir");
	system("mkdir $workdir") if (!$b);
	my $cmd = "diff ../src/dmrg $workdir/dmrg &> /dev/null";
	my $ret = system($cmd);
	system("cp -a ../src/dmrg $workdir/") if ($ret != 0);
	$cmd = "diff ../src/observe $workdir/observe &> /dev/null";
	$ret = system($cmd);
	system("cp -av  ../src/observe $workdir/") if ($ret != 0);
	chdir("$workdir/");
}
