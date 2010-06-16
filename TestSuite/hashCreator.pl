#! /usr/bin/perl -w

use strict;
use Getopt::Long;
my ($inputFile,$inputDir,$outputFile,$extension) = (undef,undef,"hashTable.txt","spec"); #Defaults for command options

#Command options for stating the file or directory to hash. Values will be stored in the output file. The extension option only applies to directories.
GetOptions("f=s" => \$inputFile, "d=s" => \$inputDir, "o=s" => \$outputFile, "x=s" => \$extension);

if(!defined($inputFile) && !defined($inputDir)) {
	print "ERROR: An input file or directory is needed!\n";
} elsif(defined($inputFile) && defined($inputDir)) {
	print "ERROR: Only one type of input can be specified!\n";
} else {
	
	if (defined($inputDir)) {
		hashDirectory($inputDir,$outputFile,$extension) if(validateDirectory($inputDir));
	} else {
		hashFile($inputFile,$outputFile) if(validateFile($inputFile));
	}
}

sub validateDirectory	#Verifies if the directory given exists
{
	my ($dir) = @_;
	
	$inputDir = "$dir/" if($dir !~ /\/$/);	#Add a backslash to the end if needed

	if(-d $inputDir) {	#Check if the directory exists
		return 1;
	} else {
		print "ERROR: The directory does not exists.\n";
		return 0;
	}
}

sub validateFile	#Verifies if the file given exists
{
	my ($file) = @_;

	if(-e $file) {		#Check if the file exists
		return 1;
	} else {
		print "ERROR: The input file does not exists.\n";
		return 0;
	}
}

sub hashDirectory	#Produce a unique hash for each file with the extension in the directory
{
	my ($dir,$output,$ext) = @_;
	my @searchFiles = glob "$dir*$ext";	#Put into an array all the files
	my $countKeys = 0;

	if (@searchFiles) {	#If at least one valid file was found
		my $file;
		my $hashKey;
		
		foreach $file(@searchFiles) {	#Create a hash value for all the files
			$countKeys += $countKeys if(hashFile($file,$output));	
		}	
		print "Created $countKeys new hash keys from ".@searchFiles." files.\n";
	} else {
		print "No files in the directory $dir with extension $ext were found.\n";
	}
}

sub hashFile	#Produce a hash value for a file
{
	my ($file,$output) = @_;
	my $hashKey;
	
	open(MYHASH,"md5sum $file|") || die "Cannot open for system command: $!\n";	#Use 'md5sum' to create 128-bit hash key
	$hashKey = substr(<MYHASH>,0,10);	#Use only the first 10 characters as the hash key
	close(MYHASH);
	
	if(!foundHashKey($hashKey,$output)) {	#If the hash value is not found in the output file
		open(MYOUTPUT,">>$output") || die "Cannot open file $output: $!\n";
		print MYOUTPUT "$hashKey\n";	#Add the hash key to the output file
		close(MYOUTPUT);
		print "A new key was created for $file\n";
		
		return 1;
	} else {
		return 0;
	}
}

sub foundHashKey	#Verify if a certain hash key is found in the output file
{
	my ($searchKey,$outfile) = @_;
	my $key;
	my $found = 0;

	open(FILE,$outfile) || die "Cannot open file $outfile: $!\n";
	while(<FILE>) {		#Iterate through the all the hash values
		chomp;
		chomp;
		$found = 1 if($_ eq $searchKey);	#Return 'true' if the hash key was matched
	}	
	close(FILE);

	return $found;
}
