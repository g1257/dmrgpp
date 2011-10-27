use Term::ANSIColor;
use strict;
=pod
run it with:
./dmrg input507-45a.inp 2>&1 | perl ~/programs/github/PsimagLite/scripts/colorOutput.pl   >& out &
tail out

reset the screen with:
perl colorOutput.pl reset

See available colors here: http://en.wikipedia.org/wiki/ANSI_escape_code

=cut

my ($needReset)=@ARGV;

if (defined($needReset)) {
	print color 'reset';
	exit(0);
}

my %colors;

$colors{"WaveFunctionTransf"}='yellow';
$colors{"Truncation"}='green';
$colors{"DmrgSolver"}='red';
$colors{"LeftRightSuper"}='magenta';
$colors{"Diag."}='cyan';
$colors{"LanczosSolver"}='bold white';
#$colors{"DensityMatrixLocal"}='darkgray';

$SIG{INT}=\&myhandler;

while(<STDIN>) {
	chomp;
	procLine($_);
	print "\n";
}

sub procLine
{
	my ($t)=@_;
	print color 'reset';
	if (!/\:/) {
		print "$t";
		return;
	}
	my @temp = split/\:/,$t;
	#print STDERR "--$temp[0]--\n";
	my $c = $colors{$temp[0]};
	if (!defined($c)) {
		print color 'reset';
		print "$t";
		return;
	}
	#print STDERR " c = $c----\n";
	print color "$c";
	print "$t";
	print color 'reset';
}



sub myhandler
{
	print color 'reset';
}

