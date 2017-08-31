#!/usr/bin/perl

use strict;
use warnings;
use Tk;
require Tk::Dialog;

my ($file,$keep,$dx,$dy) = @ARGV;
defined($file) or die "USAGE: $0 file\n";
if (defined($keep)) {
	if ($keep ne "i" && $keep ne "b") {
		die "$0: expected i (interactive) or b (batch) instead of $keep\n";
	}
} else {
	$keep = "b";
}

defined($dx) or $dx = 10;
defined($dy) or $dy = 10;

my @colors;
my ($xtotal,$ytotal,$GlobalOmegaMax) = loadColors(\@colors,$file);

my $mw = MainWindow->new;
$mw->title("Spectrum2Tk");

my @labels;

my $totalwidth=$xtotal*$dx;
my $totalheight=$ytotal*$dy;
my $xScale = 4;
my $scaleWidth = $xScale*$dx;
my $cnvMain = $mw->Canvas(-width=>$totalwidth, -height=>$totalheight);
my $cnvScale = $mw->Canvas(-width=>$scaleWidth, -height=>$totalheight);
$cnvMain->CanvasBind(-command => \&canvasClicked);

for (my $j = 0; $j < $ytotal; ++$j) {
	my $y = $j*$dy;
	for (my $i = 0; $i < $xtotal; ++$i) {
		my $x = $i*$dx;
		my $color = findColor($i,$j);
		addRectangle($x,$y,$color,$cnvMain);
	}
}

for (my $j = 0; $j < $ytotal; ++$j) {
	my $y = $j*$dy;
	my $color = rgbColor(int(($ytotal-$j-1)*255/$ytotal));
	for (my $i = 0; $i < $xScale; ++$i) {
		my $x = $i*$dx;
		addRectangle($x,$y,$color,$cnvScale);
	}
}

#menues
$mw->configure(-menu => my $mnb = $mw->Menu);
my $itemsFile = [
[qw/command ~Open -accelerator Ctrl-o -command/ => \&fileOpen],
[qw/command ~Save -accelerator Ctrl-s -command/ => \&fileSave],
'',
[qw/command ~Quit -accelerator Ctrl-q -command/ => \&exit]
];
my $mnbFile = $mnb->cascade(-label => '~File',
                            -menuitems => $itemsFile);

$cnvMain->pack(-side=>"left");
$cnvScale->pack(-side=>"right");

my $GlobalOutFile;

if ($keep eq "b") {
	$GlobalOutFile = $file;
	$GlobalOutFile =~ s/\..*$//;
	$GlobalOutFile .= ".eps";
	my $time1 = $mw->repeat(1000 => \&timer1);
}

MainLoop;

sub timer1
{
	saveToFile($GlobalOutFile,$cnvMain);
	my $scaleFilename = $GlobalOutFile;
	$scaleFilename =~ s/\..*$//;
	$scaleFilename .= "Scale.eps";
	saveToFile($scaleFilename,$cnvScale);

	$mw->withdraw();

	my $width = 7; # in cm
	my $epsFile = $GlobalOutFile;
	my $texFile = tikzCreate($epsFile,$width,$GlobalOmegaMax);
	system("pdflatex --interaction nonstopmode $texFile &> /dev/null");
	my $pdfFile = $texFile;
	$pdfFile =~ s/\..*$//;
	$pdfFile .= ".pdf";
	system("evince $pdfFile");
	exit 0;
}

sub addRectangle
{
	my ($x,$y,$color,$cnv) = @_;
	my $xp = $x + $dx;
	my $yp = $y + $dy;
	my $id = $cnv->createRectangle($x,$y,$xp,$yp,
	-fill=>$color, -width=>0,-outline =>$color);
}

sub loadColors
{
	my ($c,$file) = @_;
	open(FILE, "<", "$file") or die "$0: Cannot open file $file : $!\n";
	$_=<FILE>;
	chomp;
	my ($ytotal,$xtotal,$omegaMax) = split;
	my $counter = 0;
	while (<FILE>) {
		my @temp = split;
		$c->[$counter++] = \@temp;
	}

	close(FILE);
	print STDERR "$0: Read $counter lines, xtotal=$xtotal, ytotal=$ytotal\n";
	return ($xtotal,$ytotal,$omegaMax);
}

sub findColor
{
	my ($ind,$jnd) = @_;
	my $jnd2 = $ytotal - $jnd - 1;
	my $a = $colors[$jnd2];
	defined($a) or die "$0: findColor undefined for ind=$ind\n";
	defined($a->[$ind]) or die "$0: findColor undefined for ind=$jnd\n";
	return rgbColor($a->[$ind]);
}

sub rgbColor
{
	my ($value) = @_;
	my $green = 255-$value;
	my $blue = 255;
	my $red = $green;
	my $color = sprintf("#%02X%02X%02X", $red, $green, $blue);
	return $color;
}

sub findColor2
{
	my ($ind,$jnd) = @_;
	my $sum = $ind + $jnd;
	my $red = int(256*$ind/$xtotal);
	my $green = int(256*$sum/($xtotal+$ytotal));
	my $blue = int(256*$jnd/$ytotal);
	my $color = sprintf("#%02X%02X%02X", $red, $green, $blue);
	return $color;
}

sub fileSave
{
   my @types =
       (["PS files", [qw/.PS /]],
        ["All files",        '*'],
       );
   my $currentFile= $mw->getSaveFile(-filetypes => \@types);
   saveToFile($currentFile,$cnvMain);
}

sub saveToFile
{
	my ($currentFile,$cnv) = @_;
   $cnv->postscript(-file => $currentFile);
}

sub tikzCreate
{
	my ($graph,$w,$ymax) = @_;
	my $h = tikzHeight($w,$graph);
	my $graphScale = $graph;
	$graphScale =~ s/\..*$//;
	$graphScale .= "Scale.eps";
	my $texFile = $graph;
	$texFile =~ s/\..*$//;
	$texFile .= ".tex";
	open(FOUT, ">", "$texFile") or die "$0: Cannot write to $graph : $!\n";
	print FOUT<<EOF;
\\documentclass[border=0pt]{standalone}
\\usepackage{tikz}
\\usepackage{pgfplots}
\\pgfplotsset{width=${w}cm,compat=1.8}
\\usepgfplotslibrary{groupplots}
\\usetikzlibrary{positioning}
\\usetikzlibrary{shadows}
\\usetikzlibrary{shapes}
\\usetikzlibrary{arrows}
\\usepackage{xcolor}

\\begin{document}

\\begin{tikzpicture}
\\begin{groupplot}[group style={rows=1,columns=2}]

	\\nextgroupplot[
		axis on top,% ----
		width=${w}cm,
		height=${h}cm,
		scale only axis,
		enlargelimits=false,
		xmin=0,
		xmax=6.28,
		ymin=0,
		ymax=$ymax,
		xlabel=\$k\$,
		ylabel=\$E\$
		]

	\\addplot[thick,blue] graphics[xmin=0,ymin=0,xmax=6.28,ymax=$ymax] {$graph};

		\\nextgroupplot[
		axis on top,% ----
		width=0.5cm,
		height=${h}cm,
		scale only axis,
		enlargelimits=false,
		ymin=0,
		ymax=$ymax,
		xmin=0,
		xmax=1,
		xtickmax=-1
		]
	\\addplot[thick,blue] graphics[ymin=0,ymax=$ymax,xmin=0,xmax=1] {$graphScale};
\\end{groupplot}
\\end{tikzpicture}
\\end{document}

EOF
	close(FOUT);
	print STDERR "$0: File $texFile written\n";
	return $texFile;
}

sub tikzHeight
{
	my ($w,$graph) = @_;
	my $h;
	open(FIN, "<", $graph) or die "$0: Cannot open $graph : $!\n";
	while (<FIN>) {
		chomp;
		if (/%%BoundingBox:/) {
			my @temp = split;
			next if (scalar(@temp) != 5);
			my $xl = abs($temp[3] - $temp[1]);
			my $yl = abs($temp[4] - $temp[2]);
			$h = int($w*$yl/$xl);
		}
	}

	close(FIN);

	defined($h) or die "$0: No BoundingBox in $graph?\n";
	return $h;
}

