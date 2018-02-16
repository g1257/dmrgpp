// Created automatically by configure.pl
// DO NOT EDIT because file will be overwritten each
// time you run configure.pl with the second argument set to 1
// This file should be commited
#include "DmrgDriver1.h"


typedef PsimagLite::CrsMatrix<std::complex<RealType> > SparseMatrixInstance26Type;
typedef PsimagLite::Geometry<std::complex<RealType> ,PsimagLite::InputNg<Dmrg::InputCheck>::Readable,Dmrg::ProgramGlobals> GeometryInstance26Type;

typedef Dmrg::MatrixVectorKron<
 Dmrg::ModelBase<
  Dmrg::ModelHelperLocal<
   Dmrg::LeftRightSuper<Dmrg::BasisWithOperators<Dmrg::Operators<Dmrg::Basis<SparseMatrixInstance26Type> >  >,Dmrg::Basis<SparseMatrixInstance26Type> >
  >,
  ParametersDmrgSolverType,
  InputNgType::Readable,
  GeometryInstance26Type
 >
> MatrixVector26Type;

typedef PsimagLite::LanczosSolver<PsimagLite::ParametersForSolver<GeometryInstance26Type::RealType>,
	MatrixVector26Type, MatrixVector26Type::VectorType> LanczosSolver26Type;

template void mainLoop4<LanczosSolver26Type,Dmrg::VectorWithOffset<std::complex<RealType> > >
(LanczosSolver26Type::LanczosMatrixType::ModelType::GeometryType&,
const ParametersDmrgSolverType&,
InputNgType::Readable&,
const OperatorOptions&,
PsimagLite::String);


typedef PsimagLite::CrsMatrix<std::complex<RealType> > SparseMatrixInstance27Type;
typedef PsimagLite::Geometry<std::complex<RealType> ,PsimagLite::InputNg<Dmrg::InputCheck>::Readable,Dmrg::ProgramGlobals> GeometryInstance27Type;

typedef Dmrg::MatrixVectorOnTheFly<
 Dmrg::ModelBase<
  Dmrg::ModelHelperLocal<
   Dmrg::LeftRightSuper<Dmrg::BasisWithOperators<Dmrg::Operators<Dmrg::Basis<SparseMatrixInstance27Type> >  >,Dmrg::Basis<SparseMatrixInstance27Type> >
  >,
  ParametersDmrgSolverType,
  InputNgType::Readable,
  GeometryInstance27Type
 >
> MatrixVector27Type;

typedef PsimagLite::LanczosSolver<PsimagLite::ParametersForSolver<GeometryInstance27Type::RealType>,
	MatrixVector27Type, MatrixVector27Type::VectorType> LanczosSolver27Type;

template void mainLoop4<LanczosSolver27Type,Dmrg::VectorWithOffsets<std::complex<RealType> > >
(LanczosSolver27Type::LanczosMatrixType::ModelType::GeometryType&,
const ParametersDmrgSolverType&,
InputNgType::Readable&,
const OperatorOptions&,
PsimagLite::String);

