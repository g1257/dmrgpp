// Created automatically by configure.pl
// DO NOT EDIT because file will be overwritten each
// time you run configure.pl with the second argument set to 1
// This file should be commited
#include "DmrgDriver1.h"


typedef PsimagLite::CrsMatrix<std::complex<RealType> > SparseMatrixInstance28Type;
typedef PsimagLite::Geometry<std::complex<RealType> ,PsimagLite::InputNg<Dmrg::InputCheck>::Readable,Dmrg::ProgramGlobals> GeometryInstance28Type;

typedef Dmrg::MatrixVectorStored<
 Dmrg::ModelBase<
  Dmrg::ModelHelperLocal<
   Dmrg::LeftRightSuper<Dmrg::BasisWithOperators<Dmrg::Operators<Dmrg::Basis<SparseMatrixInstance28Type> >  >,Dmrg::Basis<SparseMatrixInstance28Type> >
  >,
  ParametersDmrgSolverType,
  InputNgType::Readable,
  GeometryInstance28Type
 >
> MatrixVector28Type;

typedef PsimagLite::LanczosSolver<PsimagLite::ParametersForSolver<GeometryInstance28Type::RealType>,
	MatrixVector28Type, MatrixVector28Type::VectorType> LanczosSolver28Type;

template void mainLoop4<LanczosSolver28Type,Dmrg::VectorWithOffsets<std::complex<RealType> > >
(LanczosSolver28Type::LanczosMatrixType::ModelType::GeometryType&,
const ParametersDmrgSolverType&,
InputNgType::Readable&,
const OperatorOptions&,
PsimagLite::String);


typedef PsimagLite::CrsMatrix<std::complex<RealType> > SparseMatrixInstance29Type;
typedef PsimagLite::Geometry<std::complex<RealType> ,PsimagLite::InputNg<Dmrg::InputCheck>::Readable,Dmrg::ProgramGlobals> GeometryInstance29Type;

typedef Dmrg::MatrixVectorKron<
 Dmrg::ModelBase<
  Dmrg::ModelHelperLocal<
   Dmrg::LeftRightSuper<Dmrg::BasisWithOperators<Dmrg::Operators<Dmrg::Basis<SparseMatrixInstance29Type> >  >,Dmrg::Basis<SparseMatrixInstance29Type> >
  >,
  ParametersDmrgSolverType,
  InputNgType::Readable,
  GeometryInstance29Type
 >
> MatrixVector29Type;

typedef PsimagLite::LanczosSolver<PsimagLite::ParametersForSolver<GeometryInstance29Type::RealType>,
	MatrixVector29Type, MatrixVector29Type::VectorType> LanczosSolver29Type;

template void mainLoop4<LanczosSolver29Type,Dmrg::VectorWithOffsets<std::complex<RealType> > >
(LanczosSolver29Type::LanczosMatrixType::ModelType::GeometryType&,
const ParametersDmrgSolverType&,
InputNgType::Readable&,
const OperatorOptions&,
PsimagLite::String);

