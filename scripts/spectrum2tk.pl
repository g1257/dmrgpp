#!/usr/bin/perl

use strict;
use warnings;
use Tk;
require Tk::Dialog;

my ($dx,$dy) = (10,10);
my ($file,$outputFile) = @ARGV;
defined($file) or die "USAGE: $0 file\n";

my @colors;
my ($xtotal,$ytotal) = loadColors(\@colors,$file);

my $mw = MainWindow->new;
$mw->title("Spectrum2Tk");

my @labels;

my $totalwidth=$xtotal*$dx;
my $totalheight=$ytotal*$dy;
my $cnv = $mw->Canvas(-width=>$totalwidth, -height=>$totalheight);
$cnv->CanvasBind(-command => \&canvasClicked);

for (my $j = 0; $j < $ytotal; ++$j) {
	my $y = $j*$dy;
	for (my $i = 0; $i < $xtotal; ++$i) {
		my $x = $i*$dx;
		my $color = findColor($i,$j);
		addRectangle($x,$y,$color);
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

$cnv->pack();

if (defined($outputFile)) {
	my $time1 = $mw->repeat(1000 => \&timer1);
}

MainLoop;

sub timer1
{
	saveToFile($outputFile);
	exit(1);
}

sub addRectangle
{
	my ($x,$y,$color) = @_;
	my $xp = $x + $dx;
	my $yp = $y + $dy;
	my $id = $cnv->createRectangle($x,$y,$xp,$yp,
	-fill=>$color,-width=>0,-outline =>$color);
}

sub loadColors
{
	my ($c,$file) = @_;
	open(FILE,"$file") or die "$0: Cannot open file $file : $!\n";
	my ($xtotal,$ytotal);
	$_=<FILE>;
	chomp;
	($ytotal,$xtotal) = split;
	my $counter = 0;
	while (<FILE>) {
		my @temp = split;
		$c->[$counter++] = \@temp;
	}

	close(FILE);
	print STDERR "$0: Read $counter lines, xtotal=$xtotal, ytotal=$ytotal\n";
	return ($xtotal,$ytotal);
}

sub findColor
{
	my ($ind,$jnd) = @_;
	my $jnd2 = $ytotal - $jnd - 1;
	my $a = $colors[$jnd2];
	defined($a) or die "$0: findColor undefined for ind=$ind\n";
	my $green = 255-$a->[$ind];
	defined($green) or die "$0: findColor undefined for ind=$jnd\n";
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
   saveToFile($currentFile);
}

sub saveToFile
{
	my ($currentFile) = @_;
   $cnv->postscript(-file => $currentFile);
}

