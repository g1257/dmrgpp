#!/usr/bin/perl6

use v6;

sub MAIN($sites, Str $type = "vertical")
{
	my %offsets = ("vertical" => 1, "horizontal" => 2, "diagonal" => 3);
	my $offset = %offsets{"$type"};
 	my Int $center = $sites div 2;

	my $x = "";
	for 0..^$sites -> Int $j {
		next if ($j +& 1);
		next if ($j + $offset >= $sites);
		$x ~= compute($center, 0, $j, 1, $offset)~",";
		$x ~= compute($center, 0, $j, 0, $offset)~",";
	}

	$x ~~ s/","$//;
	say '"' ~ $x ~ '"';
}

sub compute($ind, $s0, $jnd, $sb1, $offset)
{
	my $ss0 = ($s0) ?? "?1" !! "";
	my $ssb0 = ($s0) ?? "" !! "?1";
	my $ssb1 = ($sb1) ?? "?1" !! "";
	my $ss1 = ($sb1) ?? "" !! "?1";
	my $bind = $ind + 1;
	my $bjnd = $jnd + $offset;
	my $op1 = "c$ss0'["~$ind~"]";
	my $op2 = "c$ssb0'["~$bind~"]";
	my $op3 = "c$ssb1"~"["~$bjnd~ "]";
	my $op4 = "c$ss1"~"["~$jnd~"]";

	if ($ind < $jnd) {
		return "<gs|$op1;$op2;$op4;$op3|gs>";
	} else {
		return "<gs|$op4;$op3;$op1;$op2|gs>";
	}
}

