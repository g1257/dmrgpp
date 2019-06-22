#!/usr/bin/perl6

use v6;

my $myself = $*PROGRAM-NAME;

sub MAIN($file)
{
# For every plot
#### Read pgfplots

	my @m;
	my Int $rows = readPlot(@m, $file);

# For every plot
#### weight plot
#### write plot
	writePlot($rows, @m);
}


sub writePlot(Int $rows, @m)
{
	my Int $cols = @m.elems;

	for 0..^$cols -> Int $col {
		for 0..^$rows -> Int $row {
			my $temp = @m[$col].list[$row];
			my @temp2 = @$temp;
			@temp2.join(' ').say;
		}

		say "";
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

		if ($line eq "") {
			my $rowsFound = @thisrow.elems;
			if (defined($rows)) {
				die "$myself: column with wrong n. of rows " ~
				   " expected $rows, found $rowsFound\n"
				  if ($rows != $rowsFound);
			} else {
				$rows = $rowsFound;
			}

			my @copy = @thisrow;
			@thisrow = ();
			@m[$y++] = @copy;
			$rows = $x;
			$x = 0;
			next;
		}

		my @temp = split(/\s/, $line);
		my $total = 3;
		@temp.elems == $total or die "$myself: Wrong line $ln in $file\n";
		@thisrow[$x++] = @temp;
	}

	note "$myself: Found $rows x $y in $file\n";
	return $rows;
}



