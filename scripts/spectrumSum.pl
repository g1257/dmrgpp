#!/usr/bin/perl6

use v6;

my $self = $*PROGRAM-NAME;

sub MAIN($file1, $file2)
{
# For every plot
#### Read pgfplots

	my @m1;
	readPlot(@m1, $file1);
	my @m2;
	readPlot(@m2, $file2);

	my @msum = sumPlot(@m1, @m2);
	writePlot(@msum);
}

sub sumPlot(@m1, @m2)
{
	my Int $total = @m1.elems;
	die "$self: Plots not of equal size $total and "~@m2.elems~"\n" if ($total != @m2.elems);

	my @msum;
	for 0..^$total -> Int $ind {
		@msum[$ind] = sumThisRow(@m1[$ind], @m2[$ind]);
	}

	return @msum;
}

sub sumThisRow(@row1, @row2)
{
	my Int $total = @row1.elems;
	die "$self: Row not of equal size $total and "~@row2.elems~"\n" if ($total != @row2.elems);

	my @a;
	@a[0] = @row1[0];
	for 1..^$total -> Int $ind {
		@a[$ind] = @row1[$ind] + @row2[$ind];
	}

	return @a;
}

sub vectorEqual(@v1, @v2)
{
	my $n = @v1.elems;
	return 0 if ($n != @v2.elems);
	for 0..^$n -> Int $ind {
		next if (@v1[$ind] == @v2[$ind]);
		die "$self: @v1[$ind] different from @v2[$ind] for $ind\n";

	}

	return 1;
}
sub writePlot(@m)
{
	my Int $total = @m.elems;

	for 0..^$total -> Int $ind {
		my @thisrow = @m[$ind];
		@thisrow.join(' ').say;
	}
}

sub readPlot(@m, $file)
{
	my $input = open($file, :r);
	my @thisrow;
	my (Int $x, Int $y) = (0, 0);
	my Int $rows;

	my Int $ln = 0;
	for $file.IO.lines -> $line {
		++$ln;

		next if ($line ~~ /^ '#' /);

		my @temp = split(/\s/, $line);
		next if (@temp.elems < 2);

		@m[$y++] = @temp;
	}

	note "$self: Found $y rows in $file\n";
}



