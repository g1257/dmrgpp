#!/usr/bin/perl  -w
use strict;
use Getopt::Long;

my ($file,$noerrors,$tabLength,$fix)=("YOU MUST SPECIFY -file filename",0,4,0);
defined($file) or die "$0: Needs one argument: the filename\n";
GetOptions ("ne" => \$noerrors,    # do not display errors
	"tl=s"   => \$tabLength,      # 
	"fix"  => \$fix,
	"file=s"=>\$file);  # fix in place

my $fout = "tmp.txt";
open(FILE,$file) or die "Cannot open file $file: $!\n";
open(FOUT,">$fout") or die "Cannot open file $fout for writing: $!\n";
while(<FILE>) {
	if (/^ +/) {
		my ($spaces,$tabs)=("","");
		getSpacesAndTabs($_,\$spaces,\$tabs);
		s/(^$spaces)/$tabs/;
	}
	print FOUT;
}

close(FOUT);
close(FILE);

if ($fix) {
	system("cp $file $file.bak");
	system("cp $fout $file");
}

sub getSpacesAndTabs
{
	my ($text,$spaces,$tabs)=@_;
	$text=~/(^ +)/;
	$$spaces = $1;
	my $ll = length($$spaces);
	($ll>=$tabLength) or die "$0: File $file line $. : Found $ll spaces, expected at least $tabLength\n";
	my $t = int($ll/$tabLength);
	my $r = $ll % $tabLength;
	$ll -= $r;
	$$spaces = "";
	for (my $i=0;$i<$ll;$i++) {
		$$spaces .= " ";
	}
	$$tabs = "";
	for (my $i=0;$i<$t;$i++) {
		$$tabs .= "\t";
	}
}


		
	