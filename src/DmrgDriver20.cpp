// Created automatically by configure.pl
// DO NOT EDIT because file will be overwritten each
// time you run configure.pl with the second argument set to 1
// This file should be commited
#include "DmrgDriver1.h"


typedef PsimagLite::CrsMatrix<std::complex<RealType> > SparseMatrixInstance40Type;
typedef PsimagLite::Geometry<std::complex<RealType> ,PsimagLite::InputNg<Dmrg::InputCheck>::Readable,Dmrg::ProgramGlobals> GeometryInstance40Type;

typedef Dmrg::MatrixVectorStored<
 Dmrg::ModelBase<
  Dmrg::ModelHelperLocal<
   Dmrg::LeftRightSuper<Dmrg::BasisWithOperators<Dmrg::Operators<Dmrg::Basis<SparseMatrixInstance40Type> >  >,Dmrg::Basis<SparseMatrixInstance40Type> >
  >,
  ParametersDmrgSolverType,
  InputNgType::Readable,
  GeometryInstance40Type
 >
> MatrixVector40Type;

typedef PsimagLite::ChebyshevSolver<PsimagLite::ParametersForSolver<GeometryInstance40Type::RealType>,
	MatrixVector40Type, MatrixVector40Type::VectorType> LanczosSolver40Type;

template void mainLoop4<LanczosSolver40Type,Dmrg::VectorWithOffsets<std::complex<RealType> > >
(LanczosSolver40Type::LanczosMatrixType::ModelType::GeometryType&,
const ParametersDmrgSolverType&,
InputNgType::Readable&,
const OperatorOptions&,
PsimagLite::String);


typedef PsimagLite::CrsMatrix<std::complex<RealType> > SparseMatrixInstance41Type;
typedef PsimagLite::Geometry<std::complex<RealType> ,PsimagLite::InputNg<Dmrg::InputCheck>::Readable,Dmrg::ProgramGlobals> GeometryInstance41Type;

typedef Dmrg::MatrixVectorKron<
 Dmrg::ModelBase<
  Dmrg::ModelHelperLocal<
   Dmrg::LeftRightSuper<Dmrg::BasisWithOperators<Dmrg::Operators<Dmrg::Basis<SparseMatrixInstance41Type> >  >,Dmrg::Basis<SparseMatrixInstance41Type> >
  >,
  ParametersDmrgSolverType,
  InputNgType::Readable,
  GeometryInstance41Type
 >
> MatrixVector41Type;

typedef PsimagLite::ChebyshevSolver<PsimagLite::ParametersForSolver<GeometryInstance41Type::RealType>,
	MatrixVector41Type, MatrixVector41Type::VectorType> LanczosSolver41Type;

template void mainLoop4<LanczosSolver41Type,Dmrg::VectorWithOffsets<std::complex<RealType> > >
(LanczosSolver41Type::LanczosMatrixType::ModelType::GeometryType&,
const ParametersDmrgSolverType&,
InputNgType::Readable&,
const OperatorOptions&,
PsimagLite::String);

