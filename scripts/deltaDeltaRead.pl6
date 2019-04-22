#!/usr/bin/perl6

use v6;

my $self = $*PROGRAM-NAME;

sub MAIN(Int $sites, $file)
{
	my %h = readDeltaDelta($sites, $file);

	my @a = %h{"a"}.list;
	my @b = %h{"b"}.list;
	note "a has " ~ @a.elems ~ " and b has " ~ @b.elems;

	my Int $center = $sites div 2;

	my @value;
	my @count;
	for 0..^$sites -> Int $s {
		next if (!defined(@a[$s]));
		next if (!defined(@b[$s]));
		my $val = @a[$s] - @b[$s];
		my Int $distance = abs($s - $center) div 2;
		if (!defined(@value[$distance])) {
			@value[$distance] = $val;
			@count[$distance] = 1;
		} else {
			@value[$distance] += $val;
			++@count[$distance];
		}
	}

	my $max = @value.elems;
	for 1..^$max -> Int $distance {
		my $val = @value[$distance]/@count[$distance];
		say "$distance " ~ $val;
	}
}

sub readDeltaDelta(Int $sites, $file)
{
	my Int $center = $sites div 2;
	my $label;
	my $readNext = 0;
	my @a;
	my @b;
	for $file.IO.lines -> $line {
		next if ($line ~~ / "CmdLine" /);

		if ($readNext) {
			my @temp = split(/\s+/, $line);
			die "$self: Wrong line $line\n" if (@temp.elems != 5);
			my ($type, $site) = findTypeAndSite($label, $center);
			if ($type == 0) {
				@a[$site] = @temp[4];
			} else {
				@b[$site] = @temp[4];
			}

			$readNext = 0;
			next;
		}

		if ($line ~~ / "Fixed all sites" /) {
			$readNext = 1;
		} else {
			$readNext = 0;
		}

		if ($line ~~ /^"<gs|"/) {
			$label = $line;
			next;
		}
	}

	return { "a" => @a, "b" => @b};
}

sub findTypeAndSite($label, Int $center)
{
	#<gs|c[0];c?1[1];c'[64];c?1'[65]|gs>
	my @temp = split(/"|"/, $label);
	@temp.elems == 3 or die "$self: Wrong label $label\n";
	my @temp2 = split(/";"/, @temp[1]);
	@temp2.elems == 4 or die "$self: Wrong operators "~@temp[1]~"\n";
	my @sites;
	my @spins;
	my ($removedCenter, $removedCenterP1) = (0, 0);
	for 0..^4 -> Int $ind {
		my $op = @temp2[$ind];
		my ($spin, $site) = getSpinAndSite($op);

		if ($site == $center and !$removedCenter) {
			$removedCenter = 1;
			next;
		}

		if ($site == $center + 1 and !$removedCenterP1) {
			$removedCenterP1 = 1;
			next;
		}

		push @sites, $site;
		push @spins, $spin;
	}

	@sites.elems == 2 or die "$self: Wrong operators (sites): "~@temp[1]~"\n";
	my $s = (@sites[0] < @sites[1]) ?? @sites[0] !! @sites[1];
	my $indexForSpin0 = (@spins[0] == 0) ?? 0 !! 1;
	my $type = (@sites[$indexForSpin0] > @sites[1 - $indexForSpin0]);
	my $diff = ($type) ?? @sites[$indexForSpin0] - @sites[1-$indexForSpin0] !! @sites[1 - $indexForSpin0] - @sites[$indexForSpin0];
	#note $diff;
	return ($type, $s);
}

sub getSpinAndSite($op)
{
	my $spin = 0;
	if ($op ~~ /"?1"/) {
		$spin = 1;
	}

	my $site;
	if ($op ~~ /"[" (\d+) "]"/) {
		$site = $0;
	}

	return ($spin, $site);
}


