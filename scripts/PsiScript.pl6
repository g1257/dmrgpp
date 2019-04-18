#!/usr/bin/perl6

use v6;

my $self = $*PROGRAM-NAME;

sub MAIN(Str $name, *@args)
{
	$name ~~ / ^<[\w .]>+ $ / or die "$self: NAME $name not alphanumeric\n";
	my Str $psiPath = getPsiPath();
	my $scriptsPath = $psiPath ~ "/dmrgpp/scripts";
	my $command = findCommand($name, $scriptsPath);

	run $command, @args;
}

sub findCommand(Str $name2, Str $scriptsPath)
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

	return commandFromNameAndExtension($name, $extension, $scriptsPath);
}

sub commandFromNameAndExtension(Str $name, Str $extension, Str $scriptsPath)
{
	die "$self: scriptsPath/$name not in git index\n" unless isInGitIndex($name, $scriptsPath);

	$name.IO ~~ :r & :x or die "$self: scriptsPath/$name is not read AND exec\n";

	my $interpreter;
	$interpreter = "/usr/bin/perl6" if ($extension eq "pl6");
	$interpreter = "/usr/bin/perl" if ($extension eq "pl");
	die "$self: Cannot find interpreter\n" unless ($interpreter);

	my $argsForInterpreter = "-I $scriptsPath ";

	return "$interpreter $argsForInterpreter $scriptsPath/$name ";
}

sub isInGitIndex(Str $name, Str $scriptsPath)
{
	my $oldir = $*CWD.basename;
	chdir $scriptsPath;
	my $full = $scriptsPath ~ "/" ~ $name;
	$name ~~ / ^<[\w .]>+ $ / or die "$self: NAME $name not alphanumeric\n";
	my $capture = q:x/git ls-files $full/;
	#q:x
	return $capture.trim-trailing ~~ / "$name" /;
}

sub getPsiPath()
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

