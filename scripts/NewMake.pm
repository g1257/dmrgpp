#!/usr/bin/perl
=pod
Copyright (c) 2009-2014, UT-Battelle, LLC
All rights reserved

[PsimagLite, Version 1.0.0]

*********************************************************
THE SOFTWARE IS SUPPLIED BY THE COPYRIGHT HOLDERS AND
CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED
WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A
PARTICULAR PURPOSE ARE DISCLAIMED.

Please see full open source license included in file LICENSE.
*********************************************************

=cut
use warnings;
use strict;

package NewMake;

sub main
{
	local *FH = shift;
	my ($args, $drivers) = @_;
	die "NewMake::main(): needs drivers as 3rd argument\n" if (!defined($drivers));

	my %a = %$args;
	my $additional = $a{"additional"};
	my $additional2 = $a{"additional2"};
	my $additional3 = $a{"additional3"};
	my $additional4 = $a{"additional4"};
	my $ldflagsAdditionals = $a{"LDFLAGS"};
	my $cppflagsAdditionals = $a{"CPPFLAGS"};
	my $path = $a{"path"};
	my $code = $a{"code"};

	$additional = " " unless defined($additional);
	$additional2 = " " unless defined($additional2);
	$additional3 = "" unless defined($additional3);
	$additional4 = "" unless defined($additional4);
	$ldflagsAdditionals = "" unless defined($ldflagsAdditionals);
	$cppflagsAdditionals = "" unless defined($cppflagsAdditionals);
	$path = " " unless defined($path);

	my $allExecutables = combineAllDrivers($drivers,"");
	my $allCpps = combineAllDrivers($drivers,".cpp");

	my $configContent = getConfigContent($args->{"configFiles"}, $args->{"flavor"});

	print FH<<EOF;
# DO NOT EDIT!!! Changes will be lost. Use the PsiTag system to configure instead
# This Makefile was written by $0
# $code

$configContent

CPPFLAGS += -I$path../../PsimagLite -I$path../../PsimagLite/src -I${path}Engine
CPPFLAGS += $cppflagsAdditionals
LDFLAGS += $ldflagsAdditionals
all: $additional4 $allExecutables $additional3

EOF

	foreach my $ptr (@$drivers) {
		my $refptr = ref($ptr);
		my $oldmode = ($refptr eq "");
		my $what = ($oldmode) ? $ptr : $ptr->{"name"};
		my $aux = ($oldmode) ? 0 : $ptr->{"aux"};
		$aux = 0 if (!defined($aux));
		my $dotos = ($oldmode) ? "$what.o" : $ptr->{"dotos"};
		$dotos = "$what.o" if (!defined($dotos));

		print FH<<EOF;
$what.o: $what.cpp  Makefile $additional
	\$(CXX) \$(CPPFLAGS) -c $what.cpp

EOF

		if (!$aux) {
			# FIXME: Support many libs separated by commas here
			my $libs = ($oldmode) ? "" : $ptr->{"libs"};
			my $libs1 = "";
			my $libs2 = "";
			if (defined($libs) and $libs ne "") {
				$libs1 = "lib$libs.a";
				$libs2 = "-l$libs";
			}

			print FH<<EOF;
$what: $dotos $libs1
	\$(CXX) -o  $what $dotos \$(LDFLAGS) $libs2 \$(CPPFLAGS)
	\$(STRIP_COMMAND) $what

EOF
		}
	}

	print FH<<EOF;

$path../../PsimagLite/lib/libpsimaglite.a:
	\$(MAKE) -f Makefile -C $path../../PsimagLite/lib/

Makefile.dep: $allCpps $additional
	\$(CXX) \$(CPPFLAGS) -MM $allCpps  > Makefile.dep

clean::
	rm -f core* $allExecutables *.o *.dep $additional2

include Makefile.dep
EOF
}

sub configFilesList
{
	my ($basic, $optional) = @_;
	my @list = ($basic);
	if (defined($optional)) {
		push @list, $optional;
		print "$0: Info: Using @list\n";
		return @list;
	}

	$optional = $ENV{"PSITAG_CONFIG_USER"};
	push @list, $optional if $optional;
	print "$0: Info: Using @list\n";
	return @list;
}

sub combineAllDrivers
{
	my ($drivers,$extension) = @_;
	my $buffer = "";
	foreach my $ptr (@$drivers) {
		my $refptr = ref($ptr);
		my $oldmode = ($refptr eq "");
		my $what = ($oldmode) ? $ptr : $ptr->{"name"};
		my $aux = ($oldmode) ? 0 : $ptr->{"aux"};
		defined($aux) or $aux = 0;
		next if ($aux and $extension eq "");
		my $tmp = $what.$extension." ";
		$buffer .= $tmp;
	}

	return $buffer;
}

sub backupMakefile
{
	my ($dir) = @_;
	$dir = "." unless defined($dir);
	system("cp $dir/Makefile $dir/Makefile.bak") if (-r "$dir/Makefile");
	print STDERR "$0: Backup of $dir/Makefile in $dir/Makefile.bak\n";
}

sub getConfigContent
{
	my ($files, $flavor) = @_;
	my %tags;
	my $n = scalar(@$files);
	for (my $i = 0; $i < $n; ++$i) {
		PsiTag::main(\%tags, $files->[$i]);
	}

	flavorHelp(\%tags) if ($flavor eq "help");
	$flavor = PsiTag::canonicalTagName($tags{"default flavor"}->{"content"}) if ($flavor eq noFlavor());
	my $ptr = $tags{"flavor $flavor"};
	defined($ptr) or die "$0: No flavor named $flavor\n";
	my $configContent = $ptr->{"content"};
	defined($configContent) or die "$0: No configContent for tag $flavor\n";
	return PsiTag::unWrap(\%tags, $configContent);
}

sub noFlavor
{
	return "UNDEFINED";
}

sub flavorHelp
{
	my ($tags) = @_;
	print "$0: Available flavors are: ";
	foreach my $key (keys %$tags) {
		print "$1 " if ($key =~ /^flavor (.+)$/);
	}

	print "\n";
	exit(1);
}

1;

