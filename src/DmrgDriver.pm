#!/usr/bin/perl
use strict;
use warnings;

package DmrgDriver;

sub createTemplates
{

my $cppEach = 2;

my @lanczos = ("LanczosSolver","ChebyshevSolver");
my @matrixVector = ("MatrixVectorOnTheFly","MatrixVectorStored","MatrixVectorKron");
my @modelHelpers = ("Local","Su2");
my @vecWithOffsets = ("","s");
my @complexOrReal = ("RealType","std::complex<RealType> ");

my $cppFiles = 0;
my $counter = 0;
my $fout;
my $target = "TargetingBase";
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

					printInstance($counter,$target,$lanczos,$matrixVector,$modelHelper,
					$vecWithOffset,$complexOrNot);
					$counter++;
				}
			}
		}
	}
}

close(FOUT);

$cppFiles--;
print STDERR "$0: $counter instances and $cppFiles files\n";
return $cppFiles;
}

sub printInstance
{
	my ($counter,$target,$lanczos,$matrixVector,
	    $modelHelper,$vecWithOffset,$complexOrNot) = @_;
	my $sparseMatrix = "SparseMatrixInstance${counter}Type";
	my $realOrNotFromSparse = "${sparseMatrix}::value_type";
	my $ops = "Dmrg::Operators<Dmrg::Basis<$sparseMatrix> > ";
	my $basisWith = "Dmrg::BasisWithOperators<$ops >";
	my $basisWithout = "Dmrg::Basis<$sparseMatrix >";
	my $basis = $basisWithout;
	my $inputNg = "PsimagLite::InputNg<Dmrg::InputCheck>::Readable";
	my $geometry = "GeometryInstance${counter}Type";
	my $basisSuperBlock = "$basis";
	my $lrs = "Dmrg::LeftRightSuper<$basisWith,$basisSuperBlock >";
	my $lanczosType = "LanczosSolver${counter}Type";
	my $matrixVectorType = "MatrixVector${counter}Type";
	my $vecWithOffsetType = "Dmrg::VectorWithOffset${vecWithOffset}<$realOrNotFromSparse> ";
	print FOUT<<EOF;
#ifdef USE_COMPLEX
typedef PsimagLite::CrsMatrix<std::complex<RealType> > $sparseMatrix;
typedef PsimagLite::Geometry<$realOrNotFromSparse,$inputNg,Dmrg::ProgramGlobals> $geometry;
#else
typedef PsimagLite::CrsMatrix<$complexOrNot> $sparseMatrix;
typedef PsimagLite::Geometry<RealType,$inputNg,Dmrg::ProgramGlobals> $geometry;
#endif

typedef Dmrg::$matrixVector<
 Dmrg::ModelBase<
  Dmrg::ModelHelper$modelHelper<
   $lrs
  >,
  ParametersDmrgSolverType,
  InputNgType::Readable,
  $geometry
 >
> $matrixVectorType;

typedef PsimagLite::$lanczos<PsimagLite::ParametersForSolver<typename ${geometry}::RealType>,
	$matrixVectorType, typename ${matrixVectorType}::VectorType> $lanczosType;

template void mainLoop4<$lanczosType,$vecWithOffsetType>
(typename ${lanczosType}::LanczosMatrixType::ModelType::GeometryType&,
const ParametersDmrgSolverType&,
InputNgType::Readable&,
const OperatorOptions&,
PsimagLite::String);

EOF
}

sub isComplexOption
{
	my ($option) = @_;
	return ($option ne "RealType");
}

sub targetNotComplex
{
	my ($target) = @_;
	return ($target ne "TargetingTimeStep");
}

sub printHeader
{

	print FOUT<<EOF;
// Created automatically by $0
// DO NOT EDIT because file will be overwritten each
// time you run $0
#include "DmrgDriver1.h"

EOF
}

1;

