#!/usr/bin/perl  -w
use strict;
use Getopt::Long;

my ($file,$noerrors,$nowarn,$debug,$fix)=("YOU MUST SPECIFY -file filename",0,0,0);
defined($file) or die "$0: Needs one argument: the filename\n";
GetOptions ("ne" => \$noerrors,    # do not display errors
	"nw"   => \$nowarn,      # do not display warnings
	"fix"  => \$fix,
	"d" => \$debug,
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
my $namespacesOpen = 0;
my $lastLabel = "NO_LABEL";

my $fout = "tmp.txt";
open(FILE,$file) or die "Cannot open file $file: $!\n";
open(FOUT,">$fout") or die "Cannot open file $fout for writing: $!\n";
my $lforEcho = "";
my $appendNext = 0;
my $parensBalance = 0;

while(<FILE>) {
	if ($appendNext) {
		$lforEcho .= $_;
	} else {
		$lforEcho = $_;
	}

	if (/\\$/) {
		$appendNext = 1;
		$line++;
		next;
	} else {
		$appendNext = 0;
	}

	chomp;
	$line++;

	if ($lforEcho=~/^[\t ]+\#/) {
		die "FATAL: Line $.: preprocessor directives must start at line 0\n";
	}

	if ($lforEcho=~/^\#/) {
		print FOUT $lforEcho;
		$consecutiveEmptyLine = 0;
		next;
	}
	if (/^[\t]+DO_IF_DEBUG/) {
		print STDERR "Warning: Ignoring line $.\n" unless ($nowarn);
		print FOUT $lforEcho;
		$consecutiveEmptyLine = 0;
		next;
	}

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

	# space after if for or while or foreach 
	if (/if\(/ or /for\(/ or /while\(/ or /foreach\(/) {
		$lforEcho =~s/if\(/if \(/;
		$lforEcho =~s/for\(/for \(/;
		$lforEcho =~s/while\(/while \(/;
		$lforEcho =~s/foreach\(/foreach \(/;
		print "ERROR: Lack of space after if/for/foreach/while in line $line goes against rule -sai -saf -saw\n" unless ($noerrors);
	}

	# size_t is deprecated
	if (/size_t /) {
		print "WARNING: size_t is deprecated. Use SizeType instead\n" unless ($nowarn);
	}
	
	# std::vector, std::map and std::stack are deprecated
	my @stdclasses = ("vector","map","stack");
	foreach my $stdclass (@stdclasses) {
		if (/std::$stdclass<([^>]+)>/) {
			my $capitalclass = ucfirst($stdclass);
			print "WARNING: std::$stdclass is deprecated. Include \"$capitalclass.h\" and use [typename] PsimagLite::$capitalclass<$1>::Type instead.\n"
			unless ($nowarn);
		}
	}

	# std::string is deprecated
	my @stringclasses = ("string","istringstream","ostringstream","logic_error","range_error","runtime_error");
	foreach my $stringclass (@stringclasses) {
		if (/std::$stringclass /) {
			my $capitalclass = ucfirst($stringclass);
			if ($capitalclass =~ /_(.{1})/) {
				my $tmp = uc($1);
				$capitalclass =~ s/_(.{1})/$tmp/;
			}
			print "WARNING: std::$stringclass is deprecated. Include \"String.h\" and use PsimagLite::$capitalclass instead\n" unless ($nowarn);
		}
	}

	# line length
	$r = $_;
	$r =~ s/\t//g;
	my $ll = length($r);
	$r = $_;
	$r =~ s/[^\t]//g;
	$ll += length($r) * $sizeOfTabs;
	($ll<=90) or push @tooLongLines,($line,$ll);
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
	$o =~ s/\"[^\"]+\"//;
	$o =~ s/[^\{]//g;
	my $co = length($o);
	my $hasSemicolonAtTheEnd = 0;
	if (/\;[\t ]*$/) {
		$hasSemicolonAtTheEnd = 1;
	}

	my $parensOpen = $_;
	$parensOpen =~ s/[^\(]//g;
	$parensOpen = length($parensOpen);

	my $parensClosed = $_;
	$parensClosed =~ s/[^\)]//g;
	$parensClosed = length($parensClosed);
	$parensBalance += ($parensOpen - $parensClosed);

	my $label = getLabel($_);

	if ($co==0 and !$hasSemicolonAtTheEnd) {
		 my @whats = ("else if","if","for","foreach","else","switch","try","catch");
		 my $flagTmp = 0;
		 foreach my $what (@whats) {
			if ($label eq $what) {
				$flagTmp = 1;
				last;
			}
		}
		if ($flagTmp) {
			$parensBalance = $parensOpen - $parensClosed;
			if ($parensBalance==0) {
				$co = 1;
				$closeAfterNext++;
			}	
		}
	}

	push @mystack, $label if ($co==1);
	$c =~ s/\"[^\"]+\"//;
	$c =~ s/[^\}]//g;
	my $cc = length($c);
	if ($cc==1) {
		my $tmp = pop @mystack;
		$label = $tmp;
		if ($tmp eq "function") {
			$nextLineMustBeEmpty = 1;
		}
	}
	my $tmpLevel = $indentLevel;
	$indentLevel += ($co-$cc);
	if ($debug) {
		print STDERR "Line $.: tmplevel= $tmpLevel indentLevel= $indentLevel balance=$parensBalance co=$co cc=$cc label=$label\n";
		printList(\@mystack);
	}
	$tmpLevel = $indentLevel if ($co<$cc and !$braceAtTheEnd);
	if ($tmpLevel>0 and ($label eq "else" or $label eq "catch" or $label eq "else if") and $co>0) {
		$tmpLevel-- if ($_=~/} $label .*\{/);
	}
	
	#print "$_ ** $line ** $indentLevel \n"  if ($co<$cc and !$braceAtTheEnd);
	$braceAtTheEnd = 0;
	
	#check for public or private:
	if (/public\:/ || /private\:/ || /protected\:/ || /case [^\:]+\:/  || /default\:/) {
		$tmpLevel-- if ($tmpLevel>0);
	}

	#close namespace if needed
	if ($tmpLevel<$namespacesOpen) {
		print STDERR "Closing namespace namespacesOpen = $namespacesOpen\n";
		$namespacesOpen--; #= ($tmpLevel+1);
	}
	# check trailing whitespace
	my $w = $_;
	$w =~ /(^[ \t]*)(.*$)/;	
	checkTrailingWhite($1,$tmpLevel,$label);

	$namespacesOpen++ if ($label eq "namespace" and length($o)>0);

	print FOUT $lforEcho;
}
close(FILE);

if ($consecutiveEmptyLine!=1) {
	print "ERROR: No empty line at end of file\n" unless ($noerrors);
	print FOUT "\n";
}
close(FOUT);

printWarnings() unless ($nowarn);

if ($fix) {
	print STDERR "$0: Copying $file into $file.bak\n";
	system("cp $file $file.bak");
	print STDERR "$0: Copying $fout into $file\n";
	system("cp $fout $file");
}

sub printWarnings
{
	my $n = scalar(@tooLongLines);
	return if ($n == 0);
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
	my ($w,$inLevel,$label) = @_;
	# space before tabs is an error:
	if ($w =~ / \t/) {
		$lforEcho =~ s/(^ +)\t/\t/;
		print "ERROR: Space before tab found in line $line\n" unless ($noerrors);
	}
	my $tbs = $w;
	$tbs =~ s/[^\t]//g;
	my $ctbs = length($tbs);

	$inLevel -= $namespacesOpen;
	($ctbs==$inLevel) or die "$0: FATAL: Indentation level $inLevel tabs $ctbs, for line $line label=$label namespaces=$namespacesOpen\n";
}

sub getLabel
{
	my ($t)=@_;
	my $tt = $t;
	
	$tt =~ s/\".+\"//;
	
	my $label = $lastLabel;
	my $flag=0;
	my @whats = ("class","namespace","else if","if","for","foreach","else","enum","switch","return","try","catch");

	foreach my $what (@whats) {
		if (isOfType($tt,$what)) {
			$label = $what;
			$flag = 1;
			last;
		}
	}

	if (!$flag) {
		$label = "function";
	}

	$lastLabel = $label;
	return $label;
}

sub isOfType
{
	my ($tt,$what) = @_;
	return 1 if ($tt=~/^$what/);
	return 1 if ($tt=~/[ \t]+$what$/);
	return 1 if ($tt=~/[ \t]+$what[^a-zA-Z]/);
	return 0;
}

sub printList
{
	my ($a)= @_;
	my $n = scalar(@$a);
	for (my $i = 0; $i < $n; $i++) {
		print STDERR "$a->[$i] ";
	}
	print STDERR "\n";
}

