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
#include "DmrgDriver.h"

template<typename GeometryType,
		 typename TargettingType>
void mainLoop3(GeometryType& geometry,
			   const ParametersDmrgSolverType& dmrgSolverParams,
			   InputNgType::Readable& io,
			   const OperatorOptions& opOptions)
{
	typedef typename TargettingType::TargetParamsType TargetParamsType;
	typedef typename TargettingType::MatrixVectorType::ModelType ModelBaseType;

	//! Setup the Model
	Dmrg::ModelSelector<ModelBaseType> modelSelector(dmrgSolverParams.model);
	const ModelBaseType& model = modelSelector(dmrgSolverParams,io,geometry);

	if (opOptions.enabled) {
		operatorDriver(model,opOptions);
		return;
	}

	//! Read TimeEvolution if applicable:
	TargetParamsType tsp(io,model);

	//! Setup the dmrg solver:
	typedef Dmrg::DmrgSolver<TargettingType> SolverType;
	SolverType dmrgSolver(model,tsp,io);

	//! Calculate observables:
	dmrgSolver.main(geometry);
}

EOF
}


