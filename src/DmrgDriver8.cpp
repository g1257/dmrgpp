// Created automatically by ./configure.pl
// DO NOT EDIT because file will be overwritten each
// time you run ./configure.pl with the second argument set to 1
// This file should be commited
#include "DmrgDriver1.h"


typedef PsimagLite::CrsMatrix<std::complex<RealType> > SparseMatrixInstance16Type;
typedef Dmrg::SuperGeometry<std::complex<RealType> ,PsimagLite::InputNg<Dmrg::InputCheck>::Readable,Dmrg::ProgramGlobals> GeometryInstance16Type;

typedef Dmrg::MatrixVectorStored<
 Dmrg::ModelBase<
  Dmrg::ModelHelperLocal<
   Dmrg::LeftRightSuper<Dmrg::BasisWithOperators<Dmrg::Basis<SparseMatrixInstance16Type> ,1>,Dmrg::Basis<SparseMatrixInstance16Type> >
  >,
  ParametersDmrgSolverType,
  InputNgType::Readable,
  GeometryInstance16Type
 >
> MatrixVector16Type;

typedef PsimagLite::LanczosSolver<PsimagLite::ParametersForSolver<GeometryInstance16Type::RealType>,
	MatrixVector16Type, MatrixVector16Type::VectorType> LanczosSolver16Type;

template void mainLoop4<LanczosSolver16Type,Dmrg::VectorWithOffsets<std::complex<RealType> , Dmrg::Qn> >
(LanczosSolver16Type::MatrixType::ModelType::SuperGeometryType&,
const ParametersDmrgSolverType&,
InputNgType::Readable&,
const OperatorOptions&);


typedef PsimagLite::CrsMatrix<std::complex<RealType> > SparseMatrixInstance17Type;
typedef Dmrg::SuperGeometry<std::complex<RealType> ,PsimagLite::InputNg<Dmrg::InputCheck>::Readable,Dmrg::ProgramGlobals> GeometryInstance17Type;

typedef Dmrg::MatrixVectorKron<
 Dmrg::ModelBase<
  Dmrg::ModelHelperLocal<
   Dmrg::LeftRightSuper<Dmrg::BasisWithOperators<Dmrg::Basis<SparseMatrixInstance17Type> ,1>,Dmrg::Basis<SparseMatrixInstance17Type> >
  >,
  ParametersDmrgSolverType,
  InputNgType::Readable,
  GeometryInstance17Type
 >
> MatrixVector17Type;

typedef PsimagLite::LanczosSolver<PsimagLite::ParametersForSolver<GeometryInstance17Type::RealType>,
	MatrixVector17Type, MatrixVector17Type::VectorType> LanczosSolver17Type;

template void mainLoop4<LanczosSolver17Type,Dmrg::VectorWithOffsets<std::complex<RealType> , Dmrg::Qn> >
(LanczosSolver17Type::MatrixType::ModelType::SuperGeometryType&,
const ParametersDmrgSolverType&,
InputNgType::Readable&,
const OperatorOptions&);

