#!/usr/bin/perl

# Always use these two lines for every perl script:
use strict;
use warnings;

=pod
The line above starts a multi-line comment.
We want to load two files $file1 and $file2 containing matrices
and we want to compute the percentage difference
element by element, and then print it.
The line below ends this multi-line comment:
=cut

# This is the label we need to look for in $file1:
my $label1="OperatorC:";

# This is the label we need to look for in $file2:
my $label2="Energy=";

# The line below reads the command line in the array @ARGV:
# and places the first word in $file1 and the second in $file2
my ($file1,$file2) = @ARGV;

# These are the two matrices.
# Note that perl has no matrix data type, so we use vectors, and
# store matrix element $m($i,$j) in $matrix[$i + $j*$rows], this is called
# linear ordering
my @matrix1;
my @matrix2;


# Usually, we would use a single load function, but here
# the matrices have different formats, so we need to load functions:
my $cols1 = loadMatrix1(\@matrix1,$file1);
print STDERR "Matrix1 in file $file1 has $cols1 columns\n";

my $cols2=loadMatrix2(\@matrix2,$file2);
print STDERR "Matrix2 in file $file2 has $cols2 columns\n";

# After that we compute the percentage difference.
# But before doing this we have to get rid of the signs since
# I haven't been careful with them for cicj, sorry.
# so we construct absolute values of matrix1 into matrix1abs
my @matrix1abs;

matrixKillSign(\@matrix1abs,\@matrix1);

#and same for 2:
my @matrix2abs;
matrixKillSign(\@matrix2abs,\@matrix2);

my @matrixPd; # <-- this matrix will contain the perc. diff.

percentageDiff(\@matrixPd,\@matrix1abs,\@matrix2abs,$cols1);

# Finally we print it:
printMatrix(\@matrixPd,$cols2);

# Now we need to write the functions we used above.
# Let's start with printMatrix, which is the easiest:

sub printMatrix
{
	my ($m,$cols)=@_;  # Read the arguments
	my $col = 0; # column counter
	foreach my $x (@$m) {
		print "$x ";
		# increment column count
		$col++;
		# start a new line if new row:
		#print STDERR "$col\n";
		if ($col==$cols) {
			print "\n";
			$col=0; # and reset the column counter
		}
	}
}

sub loadMatrix1
{
	my ($m,$f)=@_; # Read the arguments
	open(FILE, "<", $f) or die "Cannot open file $f: $!\n";
	# Note that $! above contains the error message
	# Read up to let's say "OperatorC:", which is contained in $label
	while(<FILE>) {
		last if (/^$label1/);
	}
	$_=<FILE>; # read the rows and columns, like 3 6...
	#... and store them:
	my @temp=split; # this splits a space separated line into an array
	my $rows = $temp[0];
	my $cols = $temp[1];
	# now read the rest of the matrix:
	my $i = 0; # row counter
	while(<FILE>) {
		@temp=split; # this splits a space separated line into an array
		#save this row:
		my $cols = $#temp+1;
		for (my $j=0;$j<$cols;$j++) {
			$m->[$j+$i*$cols]=$temp[$j];
		}
		#increment the row
		$i++;
		# exit if we're done:
		last if ($i==$rows);
	}
	close(FILE);
	return $#temp+1; # return the number of columns
}

# Now for loading the second matrix. Again we need a different
# function because the matrix formats are different,
# compare file1.txt with file2.txt

sub loadMatrix2
{
	my ($m,$f)=@_; # Read the arguments
	open(FILE, "<", $f) or die "Cannot open file $f: $!\n";
	# Read up to let's say "Energy=", which is contained in $label2
	while(<FILE>) {
		last if (/^$label2/);
	}
	# now read the rest of the matrix:
	my $i = 0; # row counter
	my $cols=0; # number of columns
	while(<FILE>) {
		my @temp=split; # this splits a space separated line into an array
		#save this row:
		$cols = $#temp+1;
		for (my $j=0;$j<$cols;$j++) {
			$m->[$j+$i*$cols]=realPartOf($temp[$j]);
		}
		#increment the row
		$i++;
		# exit if we're done:
		last if ($i==$cols);
	}
	close(FILE);
	return $cols; # return the number of columns
}

# This little function returns $x, if given ($x,$y):
sub realPartOf
{
	my ($xy)=@_;
	$_=$xy;
	s/\(//; # kill the starting parens
	s/,.*$//; # kill everything from the comma to the end
	return $_;

}
	
sub percentageDiff
{
	my ($mpdiff,$m1,$m2,$cols1)=@_;
	# First we do mdiff=m1-m2
	my @mdiff;
	matrixDiff(\@mdiff,$m1,$m2,$cols1);
	# Then we divide mdiff by m2:
	matrixDivide($mpdiff,\@mdiff,$m2,$cols1);
	#Now we kill the sign and multiply by 100%:
	# and also adjust the precision to .xx
	my $counter=0;
	foreach my $x (@$mpdiff) {
		$mpdiff->[$counter++]=int(abs($x)*10000)/100;
	}
}

sub matrixDiff
{
	my ($mdiff,$ma,$mb,$cols)=@_;
	my $rows = $cols; # matrix is square:
	for (my $i=0;$i<$rows;$i++) {
		for (my $j=0;$j<$cols;$j++) {
			next if (!defined($ma->[$j+$i*$cols]));
			next if (!defined($mb->[$j+$i*$cols]));
			$mdiff->[$j+$i*$cols]=$ma->[$j+$i*$cols]-$mb->[$j+$i*$cols];
		}
	}			
}

sub matrixDivide
{
	my ($mdiv,$ma,$mb,$cols)=@_;
	my $rows = $cols; # matrix is square:
	for (my $i=0;$i<$rows;$i++) {
		for (my $j=0;$j<$cols;$j++) {
			next if (!defined($ma->[$j+$i*$cols]));
			next if (!defined($mb->[$j+$i*$cols]));
			next if ($mb->[$j+$i*$cols]==0);
			$mdiv->[$j+$i*$cols]=$ma->[$j+$i*$cols]/$mb->[$j+$i*$cols];
		}
	}
}

# This does: dest(i,j) = |src(i,j)| for each element i,j
sub matrixKillSign
{
	my ($dest,$src)=@_;
	my $counter=0;
	foreach my $x (@$src) {
		$dest->[$counter++]=abs($x);
	}
}



	
	
	
	
	
