#!/usr/bin/perl -w
use strict;
package AkimaSumBatch;

sub main
{
	my ($fh,$start,$end,$points,$root,$total,$ext)=@_;
	$ext = ".txt" if !defined($ext);
	#my ($start,$end,$points)=(0,3.0,100);
	my @x;
	my @s;
	my $n = 0;
	for (my $i=0;$i<$total;$i++) {
		my $file = getFile($i,$root,$end,$ext);
		my $outfile = "tmp.txt";
		print STDERR "$0: Procing $file  This is site=$i...\n";
		my $ret = system("./akimaSpline $file $start $end $points > $outfile");
		($ret ==0) or die "./akimaSpline command failed\n";
		my @y;
		$n = loadTwoColumnData(\@x,\@y,$outfile);
		sumData(\@s,\@y,$n);
	}

	printData($fh,\@x,\@s,$n);
}

sub getFile
{
	my ($i,$root,$end,$ext) = @_;
	my $j=$i;
	my $f = "$root$j$ext";
	print STDERR "Opening file $f...\n";
	open(FILE,$f) or die "$0: Cannot open file $f: $!\n";
	print STDERR "$0: in getFile: reading file $f\n";
	open(FOUT,">tmp2.txt");
	my $firstTime=1;
	my $prev = -10;
	my $val = 0;
	#while(<FILE>) {
	#	last if (/^site/);
	#}
	while(<FILE>) {
		last if (/^#/);
		
		my @temp=split;
		if ($firstTime) {
			$firstTime=0;
			if ($temp[0]!=0) {
				print FOUT "0 $temp[1]\n";
			}
		}
		# avoid duplicates
		next if ($prev>=$temp[0]);
		$prev = $temp[0];
		$val = $temp[1];
		print FOUT "$temp[0] $temp[1]\n";
	}
	if ($prev<$end) {
		print FOUT "$end $val\n";
	}
	close(FOUT);
	close(FILE);
	return "tmp2.txt";
}

sub loadTwoColumnData
{
	my ($x,$y,$file)=@_;
	my $counter=0;
	open(FILE,$file) or die "Cannot open file $file: $!\n";
	while(<FILE>) {
		next if (/^#/);
		my @temp = split;
		last if ($#temp!=1);
		$x->[$counter] = $temp[0];
		$y->[$counter] = $temp[1];
		$counter++;
	}
	close(FILE);
	return $counter;
}

sub sumData
{
	my ($s,$v,$n)=@_;
	for (my $i=0;$i<$n;$i++) {
		$s->[$i] += $v->[$i];
	}
}

sub printData
{
    my ($fh,$x,$s,$n)=@_;
	for (my $i=0;$i<$n;$i++) {
		print $fh "$x->[$i] $s->[$i]\n";
	}
}

1;

