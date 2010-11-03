#!/usr/bin/perl -w
use strict;

my ($root,$total,$ext)=@ARGV;
$ext = ".txt" if !defined($ext);
my ($start,$end,$points)=(0,5.5,100);
my @x;
my @s;
my $n = 0;
for (my $i=0;$i<$total;$i++) {
	my $file = getFile($i);
	my $outfile = "tmp.txt";
	system("./akimaSpline $file $start $end $points > $outfile");
	my @y;
	$n = loadTwoColumnData(\@x,\@y,$outfile);
	sumData(\@s,\@y,$n);
}

printData(\@x,\@s,$n);

sub getFile
{
	my $i=shift;
	my $j=$i+1;
	my $f = "$root$j$ext";
	open(FILE,$f) or die "Cannot open file $f: $!\n";
	open(FOUT,">tmp2.txt");
	my $firstTime=1;
	my $prev = -10;
	my $val = 0;
	while(<FILE>) {
		next if (/^#/);
		
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
		print FOUT;
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
        my ($x,$s,$n)=@_;
        for (my $i=0;$i<$n;$i++) {
		print "$x->[$i] $s->[$i]\n";
	}
}

