// Created automatically by configure.pl
// DO NOT EDIT because file will be overwritten each
// time you run configure.pl with the second argument set to 1
// This file should be commited
#include "DmrgDriver1.h"


typedef PsimagLite::CrsMatrix<std::complex<RealType> > SparseMatrixInstance34Type;
typedef PsimagLite::Geometry<std::complex<RealType> ,PsimagLite::InputNg<Dmrg::InputCheck>::Readable,Dmrg::ProgramGlobals> GeometryInstance34Type;

typedef Dmrg::MatrixVectorStored<
 Dmrg::ModelBase<
  Dmrg::ModelHelperSu2<
   Dmrg::LeftRightSuper<Dmrg::BasisWithOperators<Dmrg::Operators<Dmrg::Basis<SparseMatrixInstance34Type> >  >,Dmrg::Basis<SparseMatrixInstance34Type> >
  >,
  ParametersDmrgSolverType,
  InputNgType::Readable,
  GeometryInstance34Type
 >
> MatrixVector34Type;

typedef PsimagLite::LanczosSolver<PsimagLite::ParametersForSolver<GeometryInstance34Type::RealType>,
	MatrixVector34Type, MatrixVector34Type::VectorType> LanczosSolver34Type;

template void mainLoop4<LanczosSolver34Type,Dmrg::VectorWithOffsets<std::complex<RealType> , Dmrg::EffectiveQuantumNumber<RealType> > >
(LanczosSolver34Type::LanczosMatrixType::ModelType::GeometryType&,
const ParametersDmrgSolverType&,
InputNgType::Readable&,
const OperatorOptions&,
PsimagLite::String);


typedef PsimagLite::CrsMatrix<std::complex<RealType> > SparseMatrixInstance35Type;
typedef PsimagLite::Geometry<std::complex<RealType> ,PsimagLite::InputNg<Dmrg::InputCheck>::Readable,Dmrg::ProgramGlobals> GeometryInstance35Type;

typedef Dmrg::MatrixVectorKron<
 Dmrg::ModelBase<
  Dmrg::ModelHelperSu2<
   Dmrg::LeftRightSuper<Dmrg::BasisWithOperators<Dmrg::Operators<Dmrg::Basis<SparseMatrixInstance35Type> >  >,Dmrg::Basis<SparseMatrixInstance35Type> >
  >,
  ParametersDmrgSolverType,
  InputNgType::Readable,
  GeometryInstance35Type
 >
> MatrixVector35Type;

typedef PsimagLite::LanczosSolver<PsimagLite::ParametersForSolver<GeometryInstance35Type::RealType>,
	MatrixVector35Type, MatrixVector35Type::VectorType> LanczosSolver35Type;

template void mainLoop4<LanczosSolver35Type,Dmrg::VectorWithOffsets<std::complex<RealType> , Dmrg::EffectiveQuantumNumber<RealType> > >
(LanczosSolver35Type::LanczosMatrixType::ModelType::GeometryType&,
const ParametersDmrgSolverType&,
InputNgType::Readable&,
const OperatorOptions&,
PsimagLite::String);

