#!/usr/bin/perl

use strict;
use warnings;
use Getopt::Long qw(:config no_ignore_case);
use lib ".";
use Ci;

my ($valgrind,$workdir,$ranges,$regex,$su2,$info,$sOptions,$help);
my %submit;
GetOptions(
'S=s' => \$submit{"command"},
'delay=s' => \$submit{"delay"},
'valgrind=s' => \$valgrind,
'w=s' => \$workdir,
'n=s' => \$ranges,
'e=s' => \$regex,
'i=i' => \$info,
'o=s' => \$sOptions,
'su2' => \$su2,
'h' => \$help) or die "$0: Error in command line args, run with -h to display help\n";

if (defined($help)) {
	print "USAGE: $0 [options]\n";
	print "\tIf no option is given creates inputs and batches for all ";
	print "tests\n";
	print "\t-S command\n";
	print "\t\tAfter creating inputs and batches, submit them to the queue\n";
	print "\t\tusing batchDollarized.pbs as template, and\n";
	print "\t\twith command command, usually command is qsub but you can also use\n";
	print "\t\tbash to run in the command line without a batching system.\n";
	print "\t-delay delay\n";
	print "\t\tDelay in seconds between subsequent submissions.\n";
	print "\t-n n\n";
	print Ci::helpFor("-n");
	print "\t-e regex\n";
	print "Only consider input files that match regex\n";
	print Ci::helpFor("-w");
	print "\t--valgrind tool\n";
	print "\t\tRun with valgrind using tool tool\n";
	print "\t-i number\n";
	print "\t\tPrint info for test number number\n";
	print Ci::helpFor("-su2");
	print Ci::helpFor("-h");
	exit(0);
}

defined($submit{"command"}) or $submit{"command"} = "";
$submit{"PBS_O_WORKDIR"} = ($submit{"command"} eq "qsub") ? 1 : 0;
$submit{"delay"} = 5 if (!defined($submit{"delay"}));

die "$0: delay must be numeric\n" unless ($submit{"delay"} =~ /^\d+$/);

if ($submit{"delay"} < 5 and $submit{"command"} eq "qsub") {
	die "$0: Delay for qsub must be at least 5 seconds\n";
}

defined($su2) or $su2 = 0;
defined($sOptions) or $sOptions = "";

defined($valgrind) or $valgrind = "";
defined($workdir) or $workdir = "tests";

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
my @batches;

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

	my $thisInput = Ci::getInputFilename($n);
	my %keys = Ci::getInfoFromInput($thisInput, $n);
	my $isSu2 = Ci::isSu2(\%keys);
	if ($isSu2 and !$su2) {
		print STDERR "$0: WARNING: Ignored test $n ";
		print STDERR "because it's an SU(2) test and ";
		print STDERR "you did not specify -su2\n";
		next;
	}

	next unless matchesRegex($thisInput, $regex);

	my $from = getRestartFrom("$thisInput",$n);
	my @ciAnnotations = Ci::getCiAnnotations("$thisInput",$n);
	my $totalAnnotations = scalar(@ciAnnotations);

	my $whatDmrg = Ci::readAnnotationFromKey(\@ciAnnotations, "dmrg");
	my $extraCmdArgs = $sOptions."  ".findArguments($whatDmrg);
	my $cmd = getCmd($n, $valgrind, $extraCmdArgs);

	for (my $i = 0; $i < $totalAnnotations; ++$i) {
		my ($ppLabel, $w) = Ci::readAnnotationFromIndex(\@ciAnnotations, $i);
		my $x = defined($w) ? scalar(@$w) : 0;
		next if ($x == 0);
		print "|$n| has $x $ppLabel lines\n";
		next if ($ppLabel eq "dmrg");

		if ($ppLabel eq "observe") {
			$cmd .= runObserve($n, $w, $sOptions);
			next;
		}

		my $text = join ' ', map{ qq/"$_"/ } @$w;
		$text =~ s/\"\"/\"/g;
		my $pcmd = "perl ../actionsCi.pl $ppLabel $n $text\n";
		$cmd .= $pcmd;
	}

	if ($from ne "") {
		my $bi = getDeepBatchIndex($from);
		defined($batches[$bi]) or die "$0: Cannot append to empty batch\n";
		appendToBatch($batches[$bi], $cmd);
		next;
	}

	my $batch = createBatch($n, $cmd, $submit{"PBS_O_WORKDIR"});
	die "$0: Already created $batch\n" if defined($batches[$n]);
	$batches[$n] = $batch;
}

exit(0) if ($submit{"command"} eq "");

my $totalBatches = scalar(@batches);
for (my $i = 0; $i < $totalBatches; ++$i) {
	my $batch = $batches[$i];
	next if (!defined($batch));
	submitBatch(\%submit, $batch);
}

sub findArguments
{
	my ($a) = @_;
	return "" unless defined($a);
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
	my ($n, $what, $sOptions) = @_;
	my $whatN = scalar(@$what);
	my $cmd = "";
	for (my $i = 0; $i < $whatN; ++$i) {
		$cmd .= runObserveOne($n,$i,$what->[$i],$sOptions);
		$cmd .= "\n";
	}

	return $cmd;
}

sub runObserveOne
{
	my ($n, $ind, $what, $sOptions) = @_;
	# what == arguments=something
	my $args;
	if ($what =~ /^arguments=(.+$)/) {
		$args = $1;
	}

	defined($args) or die "$0: observe must have arguments\n";

	my $output = "observe$n.txt";
	unlink($output) if ($ind == 0);
	my $inputfile = Ci::getInputFilename($n);
	my $cmd = "./observe -f $inputfile $sOptions $args >> $output";
	print "|$n|: postTest $cmd\n";
	return $cmd;
}

sub matchesRegex
{
	my ($thisInput, $regex) = @_;
	defined($regex) or return 1;
	open(FILE, "<", "$thisInput") or return 1;
	my $flag = 0;
	while (<FILE>) {
		if (/$regex/) {
			$flag = 1;
			last;
		}
	}

	close(FILE);
	return $flag;
}

sub getRestartFrom
{
	my ($file,$n) = @_;
	open(FILE, "<", "$file") or die "$0: Cannot open $file : $!\n";
	my ($so, $from, $from2);
	my $restart = 0;
	while (<FILE>) {
		chomp;
		if (/SolverOptions=(.*$)/) {
			$so = $1;
			last unless ($so =~ /restart/);
			$restart = 1;
		}

		if (/RestartFilename=(.*$)/) {
			$from = $1;
		}

		if (/InfiniteLoopKeptStates=(.*$)/) {
			$from2 = $1;
		}
	}

	close(FILE);
	return "" if ($restart == 0);
	$from = $from2 unless defined($from);
	defined($from) or die "$0: restart test $n without RestartFilename\n";
	return $from;
}

sub getDeepBatchIndex
{
	my ($from) = @_;
	# $from = data??.txt or data??
	my $copy = $from;
	$copy =~ s/\.txt$//;
	$copy =~ s/^\"//;
	$copy =~ s/\";$//;
	my $n;
	if ($copy =~ /(\d+$)/) {
		$n = $1;
	}

	defined($n) or die "$0: getDeepBatchIndex $from\n";

	my $file = Ci::getInputFilename($n);
	my $from2 = getRestartFrom($file, $n);

	return ($from2 eq "") ? $n : getDeepBatchIndex($from2);
}

sub getCmd
{
	my ($n, $tool, $extraCmdArgs) = @_;
	my $valgrind = ($tool eq "") ? "" : "valgrind --tool=$tool ";
	$valgrind .= " --callgrind-out-file=callgrind$n.out " if ($tool eq "callgrind");
	my $inputfile = Ci::getInputFilename($n);
	return "$valgrind./dmrg -f $inputfile $extraCmdArgs &> output$n.txt\n\n";
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

sub appendToBatch
{
	my ($fout, $cmd) = @_;
	open(FOUT, ">>", "$fout") or die "$0: Cannot append to $fout: $!\n";
	print FOUT "$cmd\n";
	close(FOUT);

	print STDERR "$0: $fout has been written in append mode.\n";
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
	system("cp -av ../src/dmrg $workdir/");
	system("cp -av  ../src/observe $workdir/");
	chdir("$workdir/");
}


