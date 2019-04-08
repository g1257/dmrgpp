#!/usr/bin/perl6

use v6;
my $self = $*PROGRAM-NAME;

sub MAIN($file, $label, $geometry, $sites)
{
	my @m = readMatrix($file, $label);
	my @h = getDistancesAndPairs($geometry, $sites);
	my $half = @h.elems div 2;
	for 0..^$half -> Int $distance {
		
		my @pairs = @h[$distance].list;
		my $p = @pairs.elems;
		my $sum = 0;
		for 0..^$p -> Int $pind {
			my %onepair = @pairs[$pind];
			my ($i, $j) = (%onepair{"i"}, %onepair{"j"});
			#say "***$i  **** $j";
			#die "Testing\n" if ($distance == 2);
			$sum += @m[$i][$j];
		}

		say "$distance "~$sum/$p~" $p";
	}
}

sub getDistancesAndPairs($geometry, Int $sites)
{
	my @hptr;
	for 0..^$sites -> Int $i {
		for ($i+1)..^$sites -> Int $j {
			my $distance = geometryMetric($i, $j, $geometry);
			next if ($distance < 0);
			my $ij = {"i" => $i, "j" => $j};
			@hptr[$distance].append: $ij;
		}
	}

	for 0..^$sites -> Int $i {
		my $ij = {"i" => $i, "j" => $i};
		@hptr[0].append: $ij;
	}

	return @hptr;
}

sub geometryMetric($ind, $jnd, $geometry)
{
	return ($jnd - $ind) if ($geometry eq "chain");
	die "$self: Geometry $geometry not implemented\n" unless ($geometry ~~ /^"ladder"/);
	return -1 if ($ind +& 1);
	return -1 if ($jnd +& 1);
	return ($jnd - $ind) div 2;	
}


sub readMatrix($file, $label)
{
	my Int $ln = 0;
	my $fh = open $file, :r;
	while !$fh.eof {
		++$ln;
		my $line = $fh.get;
		last if ($line ~~ / "$label" /);
	}

	++$ln;
	$_ = $fh.get;
	my @temp = split(/\s/, $_);
	die "$self: Rows Cols not found, line $ln ** $_ **\n" if (@temp.elems != 2);
	my $rows = @temp[0];
	my $cols = @temp[1];

	my $row = 0;
	my @m;
	while !$fh.eof {
		++$ln;
		my $line = $fh.get;
		my @tmp = split(/\s/, $line);
		pop @tmp;
		die "$self: Wrong number of cols, line $ln, $cols not equal "~@tmp.elems~"\n" if ($cols != @tmp.elems);

		@m[$row++] = @tmp;
		last if ($row == $rows);
	}

	return @m;
}

