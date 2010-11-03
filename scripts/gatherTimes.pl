#!/usr/bin/perl -w

use strict;

my $counter = 0;
my @times;
my @values;
my $deltaT = 0.1;
while(<STDIN>) {
	next if (/^#/);
	my @temp = split;
	$times[$counter]=$temp[0];
	$values[$counter]=$temp[1];
	$counter++;
}

$counter = 0;
my @newTs = @times;
my @newVs = @values;
foreach my $t (@times)  {
	my $tx = int($t*10);
	if (($tx %2) == 1) { # it's odd then
		# does the next one exist?
		my $ind = isInVector(\@times,$#times+1,$t + $deltaT);
		if ($ind>=0) {
			# yes, then add it:
			$newVs[$ind] += $values[$counter];
			#print STDERR "Adding to $ind the value $values[$counter] for $counter, noew = $newVs[$ind]\n";
		} else {
			# no, then tranform the time:
			$newTs[$counter] = $t+ $deltaT;
			#print STDERR "Tranforming $counter into $t+ $deltaT\n";
		}
	}
	$counter++;
}

$counter = 0;
my $prevT = -1;
foreach my $t (@newTs)  {
	my $tx = int($t*10);
	if (($tx %2) == 1) {
		$counter++;
		next;
	}
	next if ($t == $prevT);
	$prevT = $t;
	print "$t $newVs[$counter++]\n";
}



sub isInVector
{
	my ($v,$n,$what)=@_;
	for (my $i=0;$i<$n;$i++) {
		return $i if ($v->[$i]==$what);
	}
	return -1;
}

