#!/usr/bin/perl
#
use strict;
use warnings;
use utf8;

my ($file, $sites, $label) = @ARGV;
defined($sites) or die "USAGE: file sites [label]\n";
defined($label) or $label = "<gs|n|gs>";

my @data;
loadData(\@data, $file);
plotData(\@data);

sub plotData
{
	my ($d) = @_;
print "\\documentclass{article}\n";
print "\\usepackage{tikz}\n";
print "\\usepackage{pgfplots}\n";
print "\\usepackage[margin=0.5cm]{geometry}\n";

print "\\begin{document}\n";
print "\\begin{tikzpicture}\n";
my $leg = 2;
my $sitesOverLeg = $sites/$leg;
my $r = 0.15;
for (my $x = 0; $x < $sitesOverLeg; ++$x) {
	my $xx = $x*0.4;
	for (my $y = 0; $y < $leg; ++$y) {
		my $id = $y + $x*2;
		my $c = int($d->[$id]*100);
		die "$0: $c > 100 or $c < 0\n" if ($c < 0 || $c > 100);
		my $yy = $y*0.5;
		print "\\draw[fill=black!$c!white,draw=none] ($xx, $yy) circle [radius=${r}];\n";
	}
}

print "%\n";

my $dx = 0.5;
my $dy = 0.5;
my $yy = 1;
for (my $i = 0; $i < 10; ++$i) {
	my $c = $i*10;
	my $xx = $i*0.5;
	print "\\draw[draw=black,fill=black!$c!white] ($xx,$yy) rectangle (";
	print "".($xx+$dx).",".($yy+$dy).");\n";
}

print <<EOF;
\\end{tikzpicture}\\vspace{1cm}

\\begin{tikzpicture}
\\begin{axis}
\\addplot+[smooth] coordinates
{
EOF
for (my $i = 0; $i < $sitesOverLeg; ++$i) {
	my $xx = $i*0.5;
	my $id = 2*$i;
	print "($xx, ".$d->[$id].") ";
}

print <<EOF;
};
\\end{axis}
\\end{tikzpicture}
\\end{document}
EOF

}

sub loadData
{
	my ($d, $file) = @_;
	open(FILE, "<", $file) or die "$0: Cannot open file $file : $!\n";
	while (<FILE>) {
		next if (/^#/);
		if (/\Q$label/) {
			chomp;
			my @temp = split;
			my $n = scalar(@temp);
			if ($n != 5) {
				print STDERR "$0: Line $.\n";
				next;
			}

			$d->[$temp[0]] = $temp[1];
		}
	}

	close(FILE);
}

