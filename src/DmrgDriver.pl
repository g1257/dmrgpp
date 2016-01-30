#!/usr/bin/perl
use strict;
use warnings;

my $cppEach = 4;

my @targets = ("TargetingGroundState","TargetingTimeStep","TargetingCorrectionVector",
"TargetingDynamic","TargetingAdaptiveDynamic","TargetingCorrection",
"MettsTargetting");
my @lanczos = ("LanczosSolver","ChebyshevSolver");
my @matrixVector = ("MatrixVectorOnTheFly","MatrixVectorStored","MatrixVectorKron");
my @modelHelpers = ("Local","Su2");
my @vecWithOffsets = ("","s");
my @complexOrReal = ("RealType","std::complex<RealType> ");

my $cppFiles = 0;
my $counter = 0;
my $fout;
foreach my $target (@targets) {
	foreach my $complexOrNot (@complexOrReal) {
		foreach my $lanczos (@lanczos) {
			foreach my $modelHelper (@modelHelpers) {
				foreach my $vecWithOffset (@vecWithOffsets) {
					foreach my $matrixVector (@matrixVector) {
						if ($counter == 0 or $counter % $cppEach == 0) {
							if ($counter > 0) {
								close(FOUT);
								print STDERR "$0: File $fout written\n";
							}

							$fout = "DmrgDriver$cppFiles.cpp";
							$cppFiles++;
							open(FOUT,"> $fout") or die "$0: Cannot write to $fout : $!\n";
							printHeader();
						}

						printInstance($target,$lanczos,$matrixVector,$modelHelper,
						$vecWithOffset,$complexOrNot);
						$counter++;
					}
				}
			}
		}
	}
}

close(FOUT);

$cppFiles--;
print STDERR "$0: $counter instances and $cppFiles files\n";

sub printInstance
{
	my ($target,$lanczos,$matrixVector,$modelHelper,$vecWithOffset,$complexOrNot) = @_;
	my $sparseMatrix = "PsimagLite::CrsMatrix<$complexOrNot> ";
	my $ops = "Dmrg::Operators<Dmrg::Basis<$sparseMatrix> > ";
	my $basisWith = "Dmrg::BasisWithOperators<$ops >";
	my $basisWithout = "Dmrg::Basis<$sparseMatrix >";
	my $basis = $basisWithout;
	my $inputNg = "PsimagLite::InputNg<Dmrg::InputCheck>::Readable";
	my $geometry = "PsimagLite::Geometry<RealType,$inputNg,Dmrg::ProgramGlobals> ";
	my $basisSuperBlock = "$basis";
	my $lrs = "Dmrg::LeftRightSuper<$basisWith,$basisSuperBlock >";
	print FOUT<<EOF;
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
   Dmrg::VectorWithOffset$vecWithOffset<$complexOrNot>
  >
 >
>
($geometry&,
const ParametersDmrgSolverType&,
InputNgType::Readable&,
const OperatorOptions&);

EOF
}

sub printHeader
{

	print FOUT<<EOF;
#include "DmrgDriver1.h"

EOF
}


