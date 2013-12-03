#!/usr/bin/perl

use strict;
use warnings;

my ($clean) = @ARGV;

if (defined($clean)) {
	die "$0: Option $clean not understood.\nUSAGE: $0 [clean]\n" unless ($clean eq "clean");
	cleanAll();
	exit(0);
}

my $paramsFile = "params.txt";

my $b = checkForFile("inputTemplate.inp",0600,"StayAlive");
if (!$b || statOfFile("inputTemplate.inp",7)==0) {
	my $someInput = getInput("Please enter the filename to use as input template: ");
	system("perl makeTemplate.pl $someInput pbsname > inputTemplate.inp");
	print STDERR "$0: WARNING file inputTemplate.inp size is 0: might want to run makeTemplate.pl manually\n"
		if (statOfFile("inputTemplate.inp",7)==0);
}

$b = checkForFile("batchTemplate.pbs",0600,"StayAlive");
if (!$b || statOfFile("batchTemplate.pbs",7)==0) {
	my $someBatch = getInput("Please enter the filename to use a PBS Batch template: ");
	my $pbsname = getInput("Please enter a root name for the #PBS -N entry: ");
	system("perl makeTemplate.pl $someBatch $pbsname > batchTemplate.pbs");
	print STDERR "$0: WARNING file batchTemplate.pbs size is 0: might want to run makeTemplate.pl manually\n"
		if (statOfFile("batchTemplate.pbs",7)==0);
}

my $useParams = "n";
my $psimagLite;

if (-r "$paramsFile") {
	print "$0: $paramsFile exits. Use it?\n";
	print "Available: y or n\n";
	print "Default: y\n";
	print "Your answer: ";
	$_ = <STDIN>;
	chomp;
	$useParams = "y" if ($_=~/^y/i or $_ eq "");
}

if ($useParams eq "n") {

	$psimagLite = getInput("Please Enter PsimagLite directory: ");
	$psimagLite =~ s/\~/$ENV{"HOME"}/;
	$psimagLite =~ s/\/$//;

	createParamsFile($paramsFile);
} else {
	$psimagLite = readLabel($paramsFile,"PsimagLite=");
}

defined($psimagLite) or die "$0: Undefined PsimagLite directory\n";

checkForFile($paramsFile,0600);	

checkForFile("$psimagLite/drivers/combineContinuedFraction",0700);
checkForFile("$psimagLite/drivers/continuedFractionCollection",0700);
checkForFile("$psimagLite/../dmrgpp/src/dmrg",0700);
checkForFile("$psimagLite/../LanczosPlusPlus/src/lanczos",0700);

systemWrapper("rsync $psimagLite/../dmrgpp/src/dmrg .");
systemWrapper("rsync $psimagLite/../LanczosPlusPlus/src/lanczos .");


sub getInput
{
	my ($txt) = @_;
	print "$txt ";
	$_=<STDIN>;
	chomp;
	return $_;
}

sub checkForFile
{
	my ($file,$perm,$stayAlive) = @_;
	
	my $mode = statOfFile($file,2);
	if (!defined($mode)) {
		die "$0: File $file does not exist\n" unless (defined($stayAlive));
		return 0;
	}

    	my $b = ($perm & $mode);
	return 1 if ($b == $perm);
	die "$0: File $file has mode $mode, but $perm expected\n" unless defined($stayAlive);
	return 0;
}


sub statOfFile
{
	#my ($dev,$ino,$mode,$nlink,$uid,$gid,$rdev,$size,
	#    $atime,$mtime,$ctime,$blksize,$blocks) = stat($file);
	my ($file,$ind) = @_;
    	my @a = stat($file);
	return $a[$ind];
}

sub systemWrapper
{
	my ($cmd) = @_;
	defined($cmd) or die "$0: systemWrapper\n";

	system($cmd);
	print "$0: Executed $cmd\n";
}

sub createParamsFile
{
	my ($file) = @_;
	my %params;
	$params{"PsimagLite"} = $psimagLite;
	my @items = ("OmegaBegin","OmegaEnd","OmegaStep","OmegaEps");
	foreach my $item (@items) {
		$params{"$item"} = getInput("Please enter a value for $item");
	}

	printHash($file,\%params);
}

sub printHash
{
	my ($file,$ptr) = @_;
	my $filebak = $file;
	$filebak =~ s/\.txt/\.bak/;
	systemWrapper("cp $file $filebak") if (-r "$file");
	my %params = %$ptr;
	open(FOUT,">$file") or die "$0: Cannot open $file\n";
	foreach my $key (keys %params) {
		print FOUT "$key=$params{$key}\n";
	}

	print FOUT "\n\n";
	close(FOUT);

	print "$0: File $file has been written.\n";
}

sub readLabel
{
	my ($file,$label) = @_;
	my $val;
	open(FILE,$file) or die "$0: Cannot open file $file : $!\n";
	while(<FILE>) {
		chomp;
		if (/^$label(.*)/) {
			$val = $1;
			last;
		}
	}

	close(FILE);

	defined($val) or die "$0: No $label in $file\n";

	return $val;
}

sub cleanAll
{
	my @files = ("lanczos","dmrg","inputTemplate.inp","batchTemplate.pbs","params.txt");
	foreach my $file (@files) {
		unlink($file);
	}
}


