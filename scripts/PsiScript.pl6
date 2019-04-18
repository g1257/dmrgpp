#!/usr/bin/perl6

use v6;

my $self = $*PROGRAM-NAME;

sub MAIN(Str $name, *@args)
{
	my Str $dmrgppPath = getPsiPath("DMRGPP");
	my $command = findCommand($name, $dmrgppPath);

	say $command ~ @args;
}

sub findCommand(Str $name2, Str $dmrgppPath)
{
	my Str $name = $name2;
	my $extension = $name.IO.extension;

	if ($extension eq "") {
		my $b6 = ($name ~ ".pl6").IO ~~ :r & :x;
		my $b5 = ($name ~ ".pl").IO ~~ :r & :x;
		if ($b5 && $b6) {
			die "$self: Failed for $name, both $b5 $b6 v6 and v5 exist\n";
		}

		if (!$b5 && !$b6) {
			die "$self: Neither $name.pl nor $name.pl6 exist\n";
		}

		$extension = "pl"  if ($b5);
		$extension = "pl6" if ($b6);
		$name ~= ".$extension" if ($extension);
	}

	die "$self: File $name with wrong extension\n" unless ($extension);

	$name.IO ~ :rx or die "$self: $name is not rx\n";
	
	my $interpreter;
	$interpreter = "/usr/bin/perl6" if ($extension eq "pl6");
	$interpreter = "/usr/bin/perl" if ($extension eq "pl");
	die "$self: Cannot find interpreter\n" unless ($interpreter);

	my $argsForInterpreter = "-I $dmrgppPath/scripts ";
	
	return "$interpreter $argsForInterpreter $dmrgppPath/scripts/$name ";
}

sub getPsiPath(Str $what)
{
	my Str $homeDir = %*ENV{"HOME"};
	my Str $file = $homeDir ~ "/.config/PsimagLite/config";
	my $psiDir;
	for $file.IO.lines {
		next if (/^\s* "#" /);
		if (/"PSI_DIR" \s* "=" \s* (.+$) /) {
			$psiDir = $0;
			last;
		}
	}

	($psiDir) or die "$self: Cannot find PSI_DIR in $file\n";
	return ($psiDir ~~ /^ "/" /) ?? $psiDir !! $homeDir ~ "/" ~ $psiDir;
}

