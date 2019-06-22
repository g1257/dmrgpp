#!/usr/bin/perl6

use v6;
my $self = $*PROGRAM-NAME;

sub MAIN($file, $label, $geometry, $sites, $label2?)
{
	my @m = readMatrix($file, $label);
	my $mode2 = "";
	my @m2;
	my @m3;
	if ($label2) {
		if ($label2 ~~ / ".cout"$ /) {
			$mode2 = "cout";
			@m3 = getDensityFromCout($label2, "<gs|n|gs>");
		} else {
			$mode2 = "diagonal";
			 @m2 = readMatrix($file, $label2);
		}
	}

	my @h = getDistancesAndPairs($geometry, $sites);
	my $half = @h.elems;
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
			$sum -= 4*@m2[$i][$i]*@m2[$j][$j] if ($mode2 eq "diagonal");
			if ($mode2 eq "cout") {
				die "$self: Undefined m3 of $i\n" if (!defined(@m3[$i]));
				die "$self: Undefined m3 of $j\n" if (!defined(@m3[$j]));
				$sum -= @m3[$i]*@m3[$j];
			}
		}

		$sum = abs($sum) if ($mode2);
		say "$distance "~$sum/$p~" $p";
	}
}

sub getDensityFromCout($file, $label)
{
	my @a;
	# 3 0.983119 0.000000 <gs|n|gs> 1.000000
	my Int $ln = 0;
	my $fh = open $file, :r;
	while !$fh.eof {
		++$ln;
		my $line = $fh.get;
		next if (/ "CmdLine" /);
		next unless ($line ~~ / "$label" /);
		my @temp = split(/\s/, $line);
		next unless (@temp.elems == 5);
		@a[@temp[0]] = @temp[1];
	}

	return @a;
}

sub getDistancesAndPairs($geometry, Int $sites)
{
	my @hptr;
	my Int $center = $sites div 2;
	for 0..^$sites -> Int $i {
		my $distance = geometryMetric($i, $center, $geometry);
		next if ($distance <= 0);
		my $imin = ($i < $center) ?? $i !! $center;
		my $j =    ($i < $center) ?? $center !! $i;
		my $ij = {"i" => $imin, "j" => $j};
		@hptr[$distance].append: $ij;
	}

	for 0..^$sites -> Int $i {
		my $ij = {"i" => $i, "j" => $i};
		@hptr[0].append: $ij;
	}

	return @hptr;
}

sub geometryMetric($ind, $jnd, $geometry)
{
	return abs($jnd - $ind) if ($geometry eq "chain");
	die "$self: Geometry $geometry not implemented\n" unless ($geometry ~~ /^"ladder"/);
	my $b1 = ($ind +& 1);
	my $b2 = ($jnd +& 1);
	return -1 if ($b1 and !$b2);
	return -1 if ($b2 and !$b1);
	return abs($jnd - $ind) div 2;
}


sub readMatrix($file, $label)
{
	my Int $ln = 0;
	my $fh = open $file, :r;
	while !$fh.eof {
		++$ln;
		my $line = $fh.get;
		next if ($line ~~ /CmdLine/);
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

