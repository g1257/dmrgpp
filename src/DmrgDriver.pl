#!/usr/bin/perl
use strict;
use warnings;

my @targets = ("TargetingGroundState","TargetingTimeStep","TargetingCorrectionVector");
my @lanczos = ("LanczosSolver","ChebyshevSolver");
my @matrixVector = ("MatrixVectorStored","MatrixVectorOnTheFly","MatrixVectorKron");
my @modelHelpers = ("Local","Su2");
my @vecWithOffsets = ("","s");

my $sparseMatrix = "PsimagLite::CrsMatrix<RealType> ";
my $ops = "Dmrg::Operators<Dmrg::Basis<$sparseMatrix> > ";
my $basisWith = "Dmrg::BasisWithOperators<$ops >";
my $basisWithout = "Dmrg::Basis<$sparseMatrix >";

printHeader();
my $counter = 0;
foreach my $target (@targets) {
	foreach my $lanczos (@lanczos) {
		foreach my $modelHelper (@modelHelpers) {
			foreach my $vecWithOffset (@vecWithOffsets) {
				foreach my $matrixVector (@matrixVector) {
					printInstance($target,$lanczos,$matrixVector,$modelHelper,$vecWithOffset);
					$counter++;
				}
			}
		}
	}
}

print STDERR "$0: $counter instances\n";

sub printInstance
{
	my ($target,$lanczos,$matrixVector,$modelHelper,$vecWithOffset) = @_;
	my $basis = $basisWithout;
	my $geometry = "PsimagLite::Geometry<RealType,InputNgType,Dmrg::ProgramGlobals> ";
	my $basisSuperBlock = "$basis";
	my $lrs = "Dmrg::LeftRightSuper<$basisWith,$basisSuperBlock >";
	print<<EOF;
template void mainLoop3<
 $geometry,
 Dmrg::$target<
  PsimagLite::$lanczos,
  Dmrg::$matrixVector<
   Dmrg::ModelBase<
    Dmrg::ModelHelper$modelHelper<
     $lrs
    >,
    ParametersDmrgSolverType,
    InputNgType::Readable,
    $geometry
   >
  >,
  Dmrg::WaveFunctionTransfFactory<
   $lrs,
   Dmrg::VectorWithOffset$vecWithOffset<RealType>
  >
 >
>
($geometry& geometry,
const ParametersDmrgSolverType& dmrgSolverParams,
InputNgType::Readable& io,
const OperatorOptions& opOptions);

EOF
}

sub printHeader
{

	print <<EOF;
#include "DmrgDriver1.h"

void usageOperator()
{
\tstd::cerr<<"USAGE is operator -f filename -F ";
\tstd::cerr<<"fermionicSign -l label [-d dof] [-s site] [-t]\\n";
}

EOF
}


