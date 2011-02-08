#!/usr/bin/perl -w
#

my ($file)=@ARGV;

my ($label1,$label2) = ("nupNdown","nUp+nDown"); 
rearrange($file,$label1);

rearrange($file,$label2);

sub rearrange
{
	my ($file,$label)=@_;
	open(FILE,$file) or die "Cannot open file $file: $!\n";

	while(<FILE>) {
		last if (/^#Using Matrix A:/);
	}
	my $x = $_;
	my $counter = 0;
	while(<FILE>) {
		$x = doOneBlock($label,$x,$counter);
		$counter++;
	}
	close(FILE);

}

sub doOneBlock 
{
	my ($label,$saved,$counter)=@_;
	while(<FILE>) {
		if (/^#/) {
			$saved .= $_;
		}
		last if (/^site /);
	}
	my $needsPrinting = 0;
	$needsPrinting = 1 if (/\Q$label/);

	print $saved if ($needsPrinting and $counter<2);

	while(<FILE>) {
		last if (/^#/);
		next if (/^Not found #FERMIONICSIGN in file/);
		next if (/^Ignore prev. error/);
		print if ($needsPrinting);
	}
	return $_;
}

