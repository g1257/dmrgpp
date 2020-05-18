#!/usr/bin/perl
# All credit goes to Nirav P.
# All errors go to G.A.

use strict;
use warnings;
use utf8;
use lib ".";
use Honeycomb;

my ($templateInput, $templateTex) = @ARGV;
defined($templateInput) or die "USAGE: $0 templateInput [templateTex]\n";

my $honey = Honeycomb::init($templateInput);

defined($honey) or die "$0: No honey\n";

createInput("test.inp", $honey->{"info"}, $templateInput);

my $params2 = {"plot" => $honey->{"plot"}};

if (defined($templateTex)) {
	createInput("test.tex", $params2, $templateTex);
} else {
	print STDERR "$0: No templateTex given, no tex created\n";
}

sub createInput
{
	my ($file, $params, $templateInput) = @_;

	open(FOUT, ">", "$file") or die "$0: Cannot write to $file\n";

	open(FILE, "<", "$templateInput") or die "$0: Cannot open $templateInput: $!\n";

	my $plot = $params->{"plot"};

	while(<FILE>) {
		next if (/^#/);
		if (/\$([a-zA-Z0-9\[\]]+)/) {
				my $name = $1;
				my $str = "\$"."params->{\"$name\"}";
				my $val = eval "$str";
				defined($val) or die "$0: Undefined substitution for $name\n";
				s/\$\Q$name/$val/g;
		}
		print FOUT;
	}

	close(FILE);
	close(FOUT);
	print STDERR "$0: File $file has been written\n";

	return $file;
}
