#!/usr/bin/perl

use strict;
use warnings;
use Getopt::Long qw(:config no_ignore_case);
use Ci;
use lib "../scripts";
use timeObservablesInSitu;
use Metts;

my ($memory,$failed,$noSu2,$help,$workdir,$golddir,$ranges,$info);
GetOptions(
'n=s' => \$ranges,
'memory=i' => \$memory,
'f' => \$failed,
'nosu2' => \$noSu2,
'i=i' => \$info,
'h' => \$help,
'g=s' => \$golddir,
'w=s' => \$workdir);

if (defined($help)) {
	print "USAGE: $0 [options]\n";
	print "\tIf no option is given examines all runs\n";
	print Ci::helpFor("-n");
	print Ci::helpFor("-w");
	print "\t-g golddir\n";
	print "\t\tUse golddir for oracles directory instead of the ";
	print "default of oldTests/\n";
	print "\t--memory mem\n";
	print "\t\tIgnore files larger than mem\n";
	print "\t-f\n";
	print "\t\tPrint info only about failed tests\n";
	print Ci::helpFor("-nosu2");
	print "\t-i number\n";
	print "\t\tPrint info for test number number\n";
	print Ci::helpFor("-h");
	exit(0);
}

defined($memory) or $memory = 50000000;
defined($failed) or $failed = 0;
defined($noSu2) or $noSu2 = 0;
defined($workdir) or $workdir = "tests";
defined($golddir) or $golddir = "oldTests";

my @tests = Ci::getTests("inputs/descriptions.txt");
my %allowedTests = Ci::getAllowedTests(\@tests);
my $total = $tests[$#tests]->{"number"};

if (defined($info)) {
	print STDERR "$0: INFO for $info\n";
	print STDERR " ".$allowedTests{$info}."\n";
	exit(0);
}

my @inRange = Ci::procRanges($ranges, $total);
my $rangesTotal = scalar(@inRange);

die "$0: No tests specified under -n\n" if ($rangesTotal == 0);

for (my $j = 0; $j < $rangesTotal; ++$j) {
        my $n = $inRange[$j];
       if (!exists($allowedTests{"$n"})) {
		#print STDERR "$0: Test $n does not exist, ignored\n";
		next;
	}

	my $isSu2 = Ci::isSu2("inputs/input$n.inp",$n);
	if ($isSu2 and $noSu2) {
		print STDERR "$0: WARNING: Ignored test $n ";
		print STDERR "because it's NOT an SU(2) test and ";
		print STDERR "you specified -nosu2\n";
		next;
        }

	procTest($n,$workdir,$golddir);

	my %ciAnnotations = Ci::getCiAnnotations("inputs/input$n.inp",$n);
	my $whatTimeInSituObs = $ciAnnotations{"getTimeObservablesInSitu"};
	my $x = scalar(@$whatTimeInSituObs);
	if ($x > 0) {
		print "|$n| has $x getTimeObservablesInSitu lines\n";
		checkTimeInSituObs($n, $whatTimeInSituObs, $workdir, $golddir);
	}

	my $metts = $ciAnnotations{"metts"};
	$x = scalar(@$metts);
	if ($x > 0) {
		print "|$n| has $x metts lines\n";
		checkMetts($n, $metts, $workdir, $golddir);
	}

	my $whatObserve = $ciAnnotations{"observe"};
	my $whatObserveN = scalar(@$whatObserve);
	if ($whatObserveN > 0) {
		print "|$n| has $whatObserveN observe lines\n";
		checkObserve($n,$workdir,$golddir);
	}

	print "-----------------------------------------------\n";
}

sub procTest
{
	my ($n,$workdir,$golddir) = @_;
	my %newValues;
	my %oldValues;
	procCout(\%newValues, $n,$workdir);
	procCout(\%oldValues, $n, $golddir);
	compareValues(\%newValues, \%oldValues, $n);
	procMemcheck($n);
}

sub procCout
{
	my ($values, $n, $dir) = @_;
	my $file = "$dir/runForinput$n.cout";
	if (!(-r "$file")) {
		print "|$n|: No $file found\n";
		return;
	}

	open(FILE, "$file") or return;
	my @energies;
	my $counter = 0;
	while (<FILE>) {
		chomp;
		if (/DMRG\+\+ version (.*$)/) {
			$values->{"version"} = $1;
			next;
		}

		if (/lowest eigenvalue= ([^ ]+) /) {
			push(@energies, $1);
			next;
		}
	}

	close(FILE);
	$values->{"energies"} = \@energies;
	$values->{"version"} = "UNDEFINED" unless (defined($values->{"version"}));

}

sub compareValues
{
	my ($newValues, $oldValues, $n) = @_;
	my $v1 = $newValues->{"version"};
	my $v2 = $oldValues->{"version"};
	defined($v1) or $v1 = "UNDEFINED";
	defined($v2) or $v2 = "UNDEFINED";
	print "|$n|: New Version $v1, Old Version $v2\n";
	my $maxEdiff = maxEnergyDiff($newValues->{"energies"}, $oldValues->{"energies"});
	print "|$n|: MaxEnergyDiff = $maxEdiff\n";
}

sub maxEnergyDiff
{
	my ($eNew, $eOld) = @_;
	return "NEW ENERGIES UNDEFINED" if (!defined($eNew));
	return "OLD ENERGIES UNDEFINED" if (!defined($eOld));
	my $n = scalar(@$eNew);
	return "ENERGY SIZES DIFFERENT" if (scalar(@$eOld) != $n);
	return "NO ENERGIES!" if ($n == 0);
	my $maxEdiff = 0;
	for (my $i = 0; $i < $n; ++$i) {
		my $tmp = abs($eNew->[$i] - $eOld->[$i]);
		$maxEdiff = $tmp if ($tmp > $maxEdiff);
 	}

	return "$maxEdiff [out of $n]";
}


sub procMemcheck
{
	my ($n) = @_;
	my $file1 = "$workdir/output$n.txt";
	my $mode = "OK";
	my $extra = "UNDEFINED";
	my %lost;
	open(FILE,$file1) or die "$0: Cannot open $file1 : $!\n";
	while (<FILE>) {
		if (/FATAL: You are requesting/ || /mandatory label/) {
			$mode = "FATAL";
			$extra = $_;
			chomp($extra);
			last;
		}

		if (/terminate called after throwing/) {
			$mode = "throw";
			while (<FILE>) {
				if (/what\(\)/) {
					$extra = $_;
					chomp($extra);
					last;
				}
			}

			last;
		}
		next unless (/^==/);
		if (/invalid/i) {
			$mode = "invalid";
			last;
		}

		if (/uninitial/i) {
			$mode = "uninitialized";
			last;
		}
		
		if (/([^ ]+) lost: ([^ ]+) bytes/) {
			my $val = $2;
			my $name = $1;
			$val =~ s/,//;
			$lost{"$name"}= $val;
			next;
		}
	}

	close(FILE);

	if ($mode eq "FATAL" || $mode eq "throw") {
		print "$0: ATTENTION TEST $n couldn't run because of $extra\n";
		return;
	}
	
	foreach my $key (keys %lost) {
		my $val = $lost{"$key"};
		next if ($val == 0);
		print "$0: ATTENTION TEST $n $key lost ".$lost{"$key"}." bytes\n";
	}

	return if ($mode eq "OK" and $failed);

	print "|$n|: output mode $mode\n";
}

sub fileSize
{
	my ($filename) = @_;
	my ($dev,$ino,$mode,$nlink,$uid,$gid,$rdev,$size,
	          $atime,$mtime,$ctime,$blksize,$blocks)
	              = stat($filename);
	return $size;
}

sub checkTimeInSituObs
{
	my ($n, $what, $workdir, $golddir) = @_;
	my $whatN = scalar(@$what);
	for (my $i = 0; $i < $whatN; ++$i) {
		my $file = "runForinput$n.cout";
		if (!(-r "$file")) {
			print STDERR "|$n|: WARNING: $file not readable\n";
		}

		my $file1 = "$workdir/timeObservablesInSitu${n}_$i.txt";
		my $file2 = "$golddir/timeObservablesInSitu${n}_$i.txt";
		print "Comparing $file1 against $file2\n";
		my %m1 = timeObservablesInSitu::load($file1);
		my %m2 = timeObservablesInSitu::load($file2);
		if (keys %m1 != 4 or keys %m2 != 4) {
			print "\tIncorrect file $file1 or $file2\n";
			return;
		}

		my $site = $m1{"site"};
		my $label = $m1{"label"};
		if ($label ne $m2{"label"}) {
			print "\tlabels NOT EQUAL\n";
			next;
		}

		if ($site ne $m2{"site"}) {
			print "\tsites NOT EQUAL\n";
			next;
		}

		if (!vectorsEqual($m1{"times"}, $m2{"times"})) {
			print "\ttimes NOT EQUAL\n";
			next;
		}

		compareMatrices($m1{"data"}, $m2{"data"});
	}
}

sub checkMetts
{
	my ($n,$what,$workdir, $golddir) = @_;
	my $whatN = scalar(@$what);
	for (my $i = 0; $i < $whatN; ++$i) {
		my $file1 = "$workdir/metts${n}_$i.txt";
		my $file2 = "$golddir/metts${n}_$i.txt";
		my %vals1 = Metts::load($file1);
		my %vals2 = Metts::load($file2);
		compareHashes(\%vals1, \%vals2);
	}
}

sub checkObserve
{
	my ($n, $workdir, $golddir) = @_;
	my $file1 = "$workdir/observe$n.txt";
	my $file2 = "$golddir/observe$n.txt";
	print "$0: Checking $file1 against $file2\n";
	my @m1 = loadObserveData($file1);
	my @m2 = loadObserveData($file2);
	compareObserveData(\@m1, \@m2);
}

sub compareObserveData
{
	my ($m1, $m2) = @_;
	my $n1 = scalar(@$m1);
	my $n2 = scalar(@$m2);
	print "work has $n1 observe matrices -- gold has $n2 observe matrices\n";
	return if ($n1 != $n2);
	for (my $i = 0; $i < $n1; ++$i) {
		compareObserveDatum($m1->[$i], $m2->[$i]);
	}
}

sub compareObserveDatum
{
	my ($h1, $h2) = @_;
	my $l1 = $h1->{"label"};
	my $l2 = $h2->{"label"};
	if ($l1 ne $l2) {
		print "Label $l1 NOT EQUAL to $l2\n";
		return;
	}

	print "Comparing $l1\n";
	compareMatrices($h1->{"data"}, $h2->{"data"});
}

sub compareMatrices
{
	my ($m1, $m2) = @_;
	if (scalar(@$m1) < 3 || scalar(@$m2) < 3) {
		print "\tMatrix TOO SMALL\n";
		return;
	}

	if ($m1->[0] != $m2->[0]) {
		print "\tRows not equal\n";
		return;
	}

	if ($m1->[1] != $m2->[1]) {
		print "\tCols not equal\n";
		return;
	}

	my $total = $m1->[0] * $m1->[1];

	my $max = 0;
	for (my $i = 0; $i < $total; ++$i) {
		my $val = abs($m1->[$i] - $m2->[$i]);
		$max = $val if ($max < $val);
	}

	print "\tMaximum difference= $max\n";
}


sub loadObserveData
{
	my ($file) = @_;
	my @m;
	my $fh;
	if (!open($fh, "$file")) {
		print "$0: File $file NOT FOUND\n";
		return @m;
	}
	
	while (<$fh>) {
		my $label = readNextLabel($fh);
		my @m1;
		my $ret = readNextMatrix($fh, \@m1);
		if ($ret ne "ok") {
			print "$0: $ret\n";
			print "$0: $label\n" if (defined($label));
			last;
		}

		my %h = ("label" => $label, "data" => \@m1);
		push @m, \%h;
	}

	close($fh);
	return @m;
}

sub readNextLabel
{
	my ($fh) = @_;
	while (<$fh>) {
		#print;
		chomp;
		last if (/^\</);
	}

	return $_;
}

sub readNextMatrix
{
	my ($fh, $m) = @_;
	$_ = <$fh>;
	defined($_) or return "eof";
	chomp;
	if (!/^\d+ /) {
		# double label, read again
		$_ = <$fh>;
		defined($_) or return "eof";
		chomp;
	}

	my @temp = split;
	
	if (scalar(@temp) != 2) {
		return "not a matrix";
	}


	my ($rows, $cols) = @temp;
	$m->[0] = $rows;
	$m->[1] = $cols;
	for (my $i = 0; $i < $rows; ++$i) {
		$_ = <$fh>;
		defined($_) or return "file ended while reading matrix";
		chomp;
		my @temp = split;
		(scalar(@temp) == $cols) or return "cols wrong for row $i";
		for (my $j = 0; $j < $cols; ++$j) {
			$m->[2 + $i +$j*$rows] = $temp[$j];
		}
	}

	return "ok";
}

sub compareHashes
{
	my ($h1, $h2) = @_;
	my $max = 0;
	foreach my $key1 (keys %$h1) {
		my $val2 = $h2->{"$key1"};
		if (!defined($val2)) {
			print "\tNot value for $key1 in hash2\n";
			next;
		}

		my $val1 = $h1->{"$key1"};
		my $d = abs($val1->[0] - $val2->[0]);
		$max = $d if ($max < $d);
	}

	print "Maximum Error $max\n";
}


			
