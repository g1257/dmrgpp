#!/usr/bin/perl  -w
use strict;

my $indentLevel = 0;
my $line = 0;
my $consecutiveEmptyLine=0;
my $closeAfterNext = 0;
my $braceAtTheEnd = 0;
my $nextLineMustBeEmpty = 0;
my @mystack;
my $sizeOfTabs = 4;
my $comment = 0;
while(<STDIN>) {
	chomp;
	$line++;
	# ignore comments
	my $cmt = "";
	$cmt = $1 if (s/([\t ]*\/\/.*$)//);
	$cmt .= $1 if s/([\t ]*\/\*.+\*\/[\t ]*)//;
	$comment = 1 if (/\/\*/);
	if (/\*\//) {
		$comment = 0;
		$cmt .= $1 if s/(^.*\*\/)//;
		$consecutiveEmptyLine = 0;
	}
	next if ($comment);
	my $r = $_;
	$r =~ s/[\t ]//g;
	if (length($r)==0) {
		# empty line:
		(length($_)==0) or print "ERROR: Empty line $line has whitespace\n";
		($consecutiveEmptyLine==0 or length($cmt)>0) or die "FATAL: Two consecutive empty lines at $line\n";
		$consecutiveEmptyLine++ if (length($cmt)==0);
		$nextLineMustBeEmpty = 0;
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
		print "ERROR: Ending whitespace in line $line\n";
	}
	
	# space after ( or before )
	if (/\( [a-zA-Z0-9]/ || /[a-zA-Z0-9] \)/) {
		print "ERROR: Space after ( or before ) in line $line goes against rule -nprs and/or -npcs\n";
	}
	# space after if for or while
	if (/if\(/ || /for\(/ or /while\(/) {
		print "ERROR: Lack of space after if/for/while in line $line goes against rule -sai -saf -saw\n";
	}
	# line length
	$r = $_;
	$r =~ s/\t//g;
	my $ll = length($r);
	$r = $_;
	$r =~ s/[^\t]//g;
	$ll += length($r) * $sizeOfTabs;
	($ll<=80) or print STDERR "WARNING: Line $line has length $ll (using tabsize=$sizeOfTabs)\n";
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
}

($consecutiveEmptyLine==1) or die "FATAL: No empty line at end of file\n";

sub checkTrailingWhite
{
	my ($w,$inLevel) = @_;
	# space before tabs is an error:
	if ($w =~ / \t/) {
		print "ERROR: Space before tab found in line $line\n";
	}
	my $tbs = $w;
	$tbs =~ s/[^\t]//g;
	my $ctbs = length($tbs);
	($ctbs==$inLevel) or die "FATAL: Indentation level $inLevel tabs $ctbs, for line $line\n";
}
