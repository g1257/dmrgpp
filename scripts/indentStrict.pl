#!/usr/bin/perl  -w
use strict;
use Getopt::Long;

my ($file,$noerrors,$nowarn,$fix)=("YOU MUST SPECIFY -file filename",0,0,0);
defined($file) or die "$0: Needs one argument: the filename\n";
GetOptions ("ne" => \$noerrors,    # do not display errors
	"nw"   => \$nowarn,      # do not display warnings
	"fix"  => \$fix,
	"file=s"=>\$file);  # fix in place

my $indentLevel = 0;
my $line = 0;
my $consecutiveEmptyLine=0;
my $closeAfterNext = 0;
my $braceAtTheEnd = 0;
my $nextLineMustBeEmpty = 0;
my @mystack;
my $sizeOfTabs = 4;
my $comment = 0;
my @tooLongLines;

my $fout = "tmp.txt";
open(FILE,$file) or die "Cannot open file $file: $!\n";
open(FOUT,">$fout") or die "Cannot open file $fout for writing: $!\n";
my $lforEcho = "";

while(<FILE>) {
	$lforEcho = $_;
	chomp;
	$line++;
	# ignore comments
	my $cmt = "";
	if (s/([\t ]*\/\/.*$)//) {
		$cmt = $1;
		$consecutiveEmptyLine = 0;
	}
	if (s/([\t ]*\/\*.+\*\/[\t ]*)//) {
		$cmt .= $1;
		$consecutiveEmptyLine = 0;
	}
	$comment = 1 if (/\/\*/);
	if (/\*\//) {
		$comment = 0;
		$cmt .= $1 if s/(^.*\*\/)//;
		$consecutiveEmptyLine = 0;
	}
	if ($comment) {
		print  FOUT $lforEcho;
		next;
	}
	my $r = $_;
	$r =~ s/[\t ]//g;
	if (length($r)==0) {
		# empty line:
		if (length($_)>0) {
			print "ERROR: Empty line $line has whitespace\n" unless ($noerrors);
			
			my $cpy = $lforEcho;
			$lforEcho = "";
			while ($cpy =~ s/(\/\*.+\*\/)//) {
				$lforEcho .= $1;
			}
			$cpy =~ s/^[\t ]*(\/\/.*$)/$1/;
			$lforEcho .= $cpy;
			$cpy =~s/[\t ]//g;
			chomp($cpy);
			if (length($cpy)==0) {
				$lforEcho ="\n";
			}
		}
		($consecutiveEmptyLine==0 or length($cmt)>0) or die "FATAL: Two consecutive empty lines at $line\n";
		$consecutiveEmptyLine++ if (length($cmt)==0);
		$nextLineMustBeEmpty = 0;
		print FOUT $lforEcho;
		next;
	} 
	if ($nextLineMustBeEmpty!=0) {
		my $r = $_;
		$r =~ s/[\t ]//g;
		($r eq "}" or $r eq "};") or die "FATAL: Line $line must be empty to comply with rule -bap\n";
	}
	$consecutiveEmptyLine=0;
	
	#check ending space
	if (/[\t ]+$/) {
		$lforEcho =~ s/[\t ]+$//;
		print "ERROR: Ending whitespace in line $line\n" unless ($noerrors);
	}
	
	# space after ( or before )
	if (/\( ([a-zA-Z0-9])/ || /([a-zA-Z0-9]) \)/) {
		$lforEcho =~s/\( ([a-zA-Z0-9])/\($1/;
		$lforEcho =~s/([a-zA-Z0-9]) \)/$1\)/;
		print "ERROR: Space after ( or before ) in line $line goes against rule -nprs and/or -npcs\n" unless ($noerrors);
	}

	# space after if for or while
	if (/if\(/ || /for\(/ or /while\(/) {
		$lforEcho =~s/if\(/if \(/;
		$lforEcho =~s/for\(/for \(/;
		$lforEcho =~s/while\(/while \(/;
		print "ERROR: Lack of space after if/for/while in line $line goes against rule -sai -saf -saw\n" unless ($noerrors);
	}
	
	# line length
	$r = $_;
	$r =~ s/\t//g;
	my $ll = length($r);
	$r = $_;
	$r =~ s/[^\t]//g;
	$ll += length($r) * $sizeOfTabs;
	($ll<=80) or push @tooLongLines,($line,$ll);
	# check indentation level
	my $o = $_;
	if ($closeAfterNext>0) {
		my $br = "";
		for (my $i=0;$i<$closeAfterNext;$i++) {
			$br .= "}";
		}
		if ($o =~ s/\;[\t ]*$/\Q$br/) {
			$braceAtTheEnd = 1;
			$closeAfterNext =0;
		}
	}
	my $c = $o;
	$o =~ s/[^\{]//g;
	my $co = length($o);
	my $hasSemicolonAtTheEnd = 0;
	if (/\;[\t ]*$/) {
		$hasSemicolonAtTheEnd = 1;
	}
	my $label = "function";
	$label = "class" if (/[^a-zA-Z]class[^a-zA-Z]/);
	$label = "if" if (/[^a-zA-Z]if[^a-zA-Z]/);
	$label = "for" if (/[^a-zA-Z]for[^a-zA-Z]/);
	$label = "else" if (/[^a-zA-Z]else[^a-zA-Z]/);
	if ($co==0 and !$hasSemicolonAtTheEnd) {
		if ($label eq "if" or $label eq "for" or $label eq "else") {
			$co = 1;
			$closeAfterNext++;
		}
	}
	push @mystack, $label if ($co==1);
	$c =~ s/[^\}]//g;
	my $cc = length($c);
	if ($cc==1) {
		my $tmp = pop @mystack;
		if ($tmp eq "function") {
			$nextLineMustBeEmpty = 1;
		}
	}
	my $tmpLevel = $indentLevel;
	$indentLevel += ($co-$cc);
	$tmpLevel = $indentLevel if ($co<$cc and !$braceAtTheEnd);
	$tmpLevel-- if ($tmpLevel>0 and $label eq "else" and $co>0);
	#print "$_ ** $line ** $indentLevel \n"  if ($co<$cc and !$braceAtTheEnd);
	$braceAtTheEnd = 0;
	
	#check for public or private:
	if (/public\:/ || /private\:/ || /protected\:/) {
		$tmpLevel-- if ($tmpLevel>0);
	}
	# check trailing whitespace
	my $w = $_;
	$w =~ /(^[ \t]*)(.*$)/;	
	checkTrailingWhite($1,$tmpLevel);
	print FOUT $lforEcho;
}
close(FILE);

if ($consecutiveEmptyLine!=1) {
	print "ERROR: No empty line at end of file\n" unless ($noerrors);
	print FOUT "\n";
}
close(FOUT);

printWarnings() unless ($nowarn);

sub printWarnings
{
	my $c = 0;
	print "LINE\tLENGTH (tabsize=$sizeOfTabs)\n";
	foreach (@tooLongLines) {
		if ($c%2==0) {
			print "$_ ";
		} else {
			print "$_\n";
		}
		$c++;
	}
}

sub checkTrailingWhite
{
	my ($w,$inLevel) = @_;
	# space before tabs is an error:
	if ($w =~ / \t/) {
		$lforEcho =~ s/(^ +)\t/\t/;
		print "ERROR: Space before tab found in line $line\n" unless ($noerrors);
	}
	my $tbs = $w;
	$tbs =~ s/[^\t]//g;
	my $ctbs = length($tbs);
	($ctbs==$inLevel) or die "FATAL: Indentation level $inLevel tabs $ctbs, for line $line\n";
}
