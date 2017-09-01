#!/usr/bin/perl

use strict;
use warnings;
use Getopt::Long qw(:config no_ignore_case);
use Ci;
use lib "../scripts";
use timeObservablesInSitu;
use CollectBrakets;
use Metts;

my ($valgrind,$workdir,$restart,$ranges,$postprocess,$noSu2,$info,$help);
my %submit;
GetOptions(
'S=s' => \$submit{"command"},
'delay=s' => \$submit{"delay"},
'valgrind=s' => \$valgrind,
'w=s' => \$workdir,
'R' => \$restart,
'P' => \$postprocess,
'n=s' => \$ranges,
'i=i' => \$info,
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
	print "\t-i number\n";
	print "\t\tPrint info for test number number\n";
	print Ci::helpFor("-nosu2");
	print Ci::helpFor("-h");
	exit(0);
}

defined($submit{"command"}) or $submit{"command"} = "";
$submit{"PBS_O_WORKDIR"} = ($submit{"command"} eq "qsub") ? 1 : 0;

defined($noSu2) or $noSu2 = 0;

defined($valgrind) or $valgrind = "";
defined($workdir) or $workdir = "tests";
defined($restart) or $restart = 0;
defined($postprocess) or $postprocess = 0;

my $templateBatch = "batchDollarized.pbs";

prepareDir();

my @tests = Ci::getTests("../inputs/descriptions.txt");
my %allowedTests = Ci::getAllowedTests(\@tests);
my $total = $tests[$#tests]->{"number"};

if (defined($info)) {
	my $desc = $allowedTests{$info};
	defined($desc) or die "$0: No test $info\n";
	print STDERR "$0: INFO for $info\n";
	print STDERR " ".$desc."\n";
	exit(0);
}

my @inRange = Ci::procRanges($ranges, $total);
my $rangesTotal = scalar(@inRange);

die "$0: No tests specified under -n\n" if ($rangesTotal == 0);
#print STDERR "@inRange"."\n";
my $nonExistent = "";

for (my $j = 0; $j < $rangesTotal; ++$j) {
	my $n = $inRange[$j];
	if (!exists($allowedTests{"$n"})) {
		$nonExistent .= "$n ";
		next;
	}

	if ($nonExistent ne "") {
		$nonExistent = Ci::compactList($nonExistent);
		#print STDERR "$0: Test(s) $nonExistent do(es) not exist, ignored\n";
		$nonExistent = "";
	}

	my $ir = isRestart("../inputs/input$n.inp",$n);
	if ($restart and !$ir) {
		print STDERR "$0: WARNING: Ignored test $n ";
		print STDERR "because it's NOT a restart test and ";
		print STDERR "you specified -R\n";
		next;
	}

	if (!$restart and $ir) {
		print STDERR "$0: WARNING: Ignored test $n ";
		print STDERR "because it's a restart test and ";
		print STDERR "you NOT specify -R\n";
		next;
	}

	my $isSu2 = Ci::isSu2("../inputs/input$n.inp",$n);
	if ($isSu2 and $noSu2) {
		print STDERR "$0: WARNING: Ignored test $n ";
		print STDERR "because it's NOT an SU(2) test and ";
		print STDERR "you specified -nosu2\n";
		next;
	}

	my %ciAnnotations = Ci::getCiAnnotations("../inputs/input$n.inp",$n);
	my $whatTimeInSituObs = $ciAnnotations{"getTimeObservablesInSitu"};
	my $x = scalar(@$whatTimeInSituObs);
	if ($x > 0) {
		print "|$n| has $x getTimeObservablesInSitu lines\n";
		if ($postprocess) {
			runTimeInSituObs($n, $whatTimeInSituObs, \%submit);
		}
	}

	my $brakets = $ciAnnotations{"CollectBrakets"};
	$x = scalar(@$brakets);
	if ($x > 0) {
		print "|$n| has $x CollectBrakets lines\n";
		if ($postprocess) {
			runCollectBrakets($n, $brakets, \%submit);
		}
	}

	my $metts = $ciAnnotations{"metts"};
	$x = scalar(@$metts);
	if ($x > 0) {
		print "|$n| has $x metts lines\n";
		if ($postprocess) {
			runMetts($n, $metts, \%submit);
		}
	}

	my $whatObserve = $ciAnnotations{"observe"};
	my $whatObserveN = scalar(@$whatObserve);
	if ($whatObserveN > 0) {
		print "|$n| has $whatObserveN observe lines\n";
		if ($postprocess) {
			runObserve($n,$whatObserve,\%submit);
		}
	}

	next if ($postprocess);

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

sub runTimeInSituObs
{
	my ($n,$what,$submit) = @_;
	my $whatN = scalar(@$what);
	for (my $i = 0; $i < $whatN; ++$i) {
		my $file = "runForinput$n.cout";
		if (!(-r "$file")) {
			print STDERR "|$n|: WARNING: $file not readable\n";
			next;
		}

		my @temp = split(/ /, $what->[$i]);
		(scalar(@temp) == 2) or die "$0: FATAL annotation ".$what->[$i]."\n";
		my ($site, $label) = @temp;		
		my $fin;
		open($fin, "<", $file) or die "$0: Could not open $file : $!\n";
		my $fout;
		my $foutname = "timeObservablesInSitu${n}_$i.txt";
		if (!open($fout, ">", "$foutname")) {
			close($fin);
			die "$0: Could not write to $foutname: $!\n";
		}

		if ($submit->{"command"} ne "") {
			timeObservablesInSitu::main($site, $label, $fin, $fout);
		} else {
			print STDERR "|$n|: Dry run $site $label\n";
		}

		close($fin);
		close($fout);
	}
}

sub runCollectBrakets
{
	my ($n,$what,$submit) = @_;
	my $whatN = scalar(@$what);
	for (my $i = 0; $i < $whatN; ++$i) {
		my $file = "runForinput$n.cout";
		if (!(-r "$file")) {
			print STDERR "|$n|: WARNING: $file not readable\n";
			next;
		}

		#my @temp = split(/ /, $what->[$i]); # arguments ot CollectBrakets
		my $foutname = "CollectBrakets${n}_$i.txt";
		
		if ($submit->{"command"} ne "") {
			CollectBrakets::main($file, $foutname);
		} else {
			print STDERR "|$n|: Dry run $file $foutname\n";
		}
	}
}

sub runMetts
{
	my ($n,$what,$submit) = @_;
	my $whatN = scalar(@$what);
	my %actions = ("Energy" => \&Metts::energy,
	               "Density" => \&Metts::density);
	for (my $i = 0; $i < $whatN; ++$i) {
		my $file = "runForinput$n.cout";
		if (!(-r "$file")) {
			print STDERR "|$n|: WARNING: $file not readable\n";
			next;
		}

		my @temp = split(/ /, $what->[$i]);
		(scalar(@temp) == 3) or next;

		my ($label, $arg0, $arg1) = @temp;
		if (($label ne "Energy") and ($label ne "Density")) {
			die "$0: Wrong annotation: $what->[$i]\n";
		}

		my $fin;
		open($fin, "<", $file) or die "$0: Could not open $file : $!\n";
		my $fout;
		my $foutname = "metts${n}_$i.txt";
		if (!open($fout, ">", "$foutname")) {
			close($fin);
			die "$0: Could not write to $foutname: $!\n";
		}

		if ($submit->{"command"} ne "") {
			my ($sum, $counter) = $actions{"$label"}($arg0, $arg1,0, $fin);
			print $fout "#"."$label=$sum $counter\n";
		} else {
			print STDERR "|$n|: Dry run $arg0 $arg1\n";
		}

		close($fin);
		close($fout);
	}
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

	my $batch = createBatch($n,$cmd, $submit->{"PBS_O_WORKDIR"});
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

sub isRestart
{
	my ($file,$n) = @_;
	open(FILE, "<", "$file") or return 0;
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
	my $batch = createBatch($n, $cmd, $submit->{"PBS_O_WORKDIR"});
	submitBatch($submit, $batch) if ($submit->{"command"} ne "");
}

sub createBatch
{
	my ($ind,$cmd,$pbsOworkDir) = @_;
	my $file = "Batch$ind.pbs";
	open(FOUT, ">", "$file") or die "$0: Cannot write to $file: $!\n";

	open(FILE, "<", "../$templateBatch") or die "$0: Cannot open ../$templateBatch: $!\n";

	while (<FILE>) {
		s/\$PBS_O_WORKDIR/\./g if (!$pbsOworkDir);
		while (/\$\$([a-zA-Z0-9\[\]]+)/) {
			my $line = $_;
			my $name = $1;
			my $str = "\$".$name;
			my $val = eval "$str";
			defined($val) or die "$0: Undefined substitution for $name\n";
			$line =~ s/\$\$$name/$val/;
			$_ = $line;
		}

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
