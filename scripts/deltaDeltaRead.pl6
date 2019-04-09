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
	for 0..^$sites -> Int $s {
		next if (!defined(@a[$s]));
		next if (!defined(@b[$s]));
		my $val = @a[$s] - @b[$s];
		my Int $distance = abs($s - $center) div 2;
		if (defined(@value[$distance])) {
			@value[$distance] = $val;
		} else {
			@value[$distance] += $val;
		}
	}

	my $max = @value.elems;
	for 1..^$max -> Int $distance {
		say "$distance " ~ abs(@value[$distance]);
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
	my $type;
	my @sites;
	my $typeCenter;
	for 0..^4 -> Int $ind {
		my $op = @temp2[$ind];
		my ($spin, $site) = getSpinAndSite($op);
		$typeCenter = ($site +& 1) if ($spin == 0);
		next if ($site == $center or $site == $center + 1);
		push @sites, $site;
		next if ($spin != 0);
		$type = ($site +& 1);
	}

	@sites = ($center, $center + 1) if (@sites.elems == 0);
	@sites.elems == 2 or die "$self: Wrong operators (sites): "~@temp[1]~"\n";
	my $s = (@sites[0] < @sites[1]) ?? @sites[0] !! @sites[1];
	defined($type) or $type = $typeCenter;
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


