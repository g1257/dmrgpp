#!/usr/bin/perl
use strict;
use warnings;

package DmrgDriver;

sub createTemplates
{
	my ($generateSources) = @_;

my $cppEach = 2;

my @lanczos = ("LanczosSolver","ChebyshevSolver");
my @matrixVector = ("MatrixVectorOnTheFly","MatrixVectorStored","MatrixVectorKron");
my @vecWithOffsets = ("","s");
my @complexOrReal = (0, 3); 

my @values;
my $cppFiles = 0;
my $counter = 0;
my $fout;
my $target = "TargetingBase";
my @su2files;
foreach my $complexOrNot (@complexOrReal) {
	foreach my $lanczos (@lanczos) {
		foreach my $vecWithOffset (@vecWithOffsets) {
			foreach my $matrixVector (@matrixVector) {
				if ($counter == 0 or $counter % $cppEach == 0) {
					if ($counter > 0 and $generateSources) {
						close(FOUT);
						print STDERR "$0: File $fout written\n";
					}

					$fout = "DmrgDriver$cppFiles.cpp";
					$su2files[$cppFiles++] = 0;

					if ($generateSources) {
						open(FOUT,"> $fout") or die "$0: Cannot write to $fout : $!\n";
						printHeader();
					}
				}

				if ($generateSources) {
					printInstance($counter,$target,$lanczos,$matrixVector,
				$vecWithOffset,$complexOrNot,\@values);
				}

				$counter++;
			}
		}
	}
}

if ($generateSources) {
	close(FOUT);

	print STDERR "$0: $counter instances and $cppFiles files\n";
}

return @su2files;
}

sub printInstance
{
	my ($counter,$target,$lanczos,$matrixVector,$vecWithOffset,$complexOrNot,
	$values) = @_;
	my $sparseMatrix = "SparseMatrixInstance${counter}Type";
	my $basisWithout = "Dmrg::Basis<$sparseMatrix>";
	my $basisWith = "Dmrg::BasisWithOperators<$basisWithout >";
	my $basis = $basisWithout;
	my $inputNg = "PsimagLite::InputNg<Dmrg::InputCheck>::Readable";
	my $geometry = "GeometryInstance${counter}Type";
	my $basisSuperBlock = "$basis";
	my $lrs = "Dmrg::LeftRightSuper<$basisWith,$basisSuperBlock >";
	my $lanczosType = "LanczosSolver${counter}Type";
	my $matrixVectorType = "MatrixVector${counter}Type";
	my ($complexOrReal1, $complexOrReal2) = getDualRealComplex($complexOrNot);
	my $vecWithOffsetType = "Dmrg::VectorWithOffset${vecWithOffset}<$complexOrReal2, Dmrg::Qn> ";

	print FOUT<<EOF;

typedef PsimagLite::CrsMatrix<$complexOrReal1> $sparseMatrix;
typedef Dmrg::SuperGeometry<$complexOrReal1,$inputNg,Dmrg::ProgramGlobals> $geometry;

typedef Dmrg::$matrixVector<
 Dmrg::ModelBase<
  Dmrg::ModelHelperLocal<
   $lrs
  >,
  ParametersDmrgSolverType,
  InputNgType::Readable,
  $geometry
 >
> $matrixVectorType;

typedef PsimagLite::$lanczos<PsimagLite::ParametersForSolver<${geometry}::RealType>,
	$matrixVectorType, ${matrixVectorType}::VectorType> $lanczosType;

template void mainLoop4<$lanczosType,$vecWithOffsetType>
(${lanczosType}::MatrixType::ModelType::SuperGeometryType&,
const ParametersDmrgSolverType&,
InputNgType::Readable&,
const OperatorOptions&);

EOF

my $value = "Local$complexOrNot";

my $seen  = 0;
foreach my $item (@$values) {
	if ($item eq $value) {
		$seen = 1;
		last;
	}
}

return if ($seen);

push @$values, $value;
foreach my $item (@$values) {
	system("echo \'$item\' >> out.txt");
}

system("echo \"-------------------\" >> out.txt");

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
// time you run $0 with the second argument set to 1
// This file should be commited
#include "DmrgDriver1.h"

EOF
}

sub getDualRealComplex
{
	my ($n) = @_;
	return (getOneRealComplex($n & 1), getOneRealComplex($n & 2));
}

sub getOneRealComplex
{
	my ($n) = @_;
	return ($n == 0) ? "RealType" : "std::complex<RealType> ";
}

1;

