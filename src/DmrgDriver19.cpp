// Created automatically by configure.pl
// DO NOT EDIT because file will be overwritten each
// time you run configure.pl with the second argument set to 1
// This file should be commited
#include "DmrgDriver1.h"


typedef PsimagLite::CrsMatrix<std::complex<RealType> > SparseMatrixInstance38Type;
typedef PsimagLite::Geometry<std::complex<RealType> ,PsimagLite::InputNg<Dmrg::InputCheck>::Readable,Dmrg::ProgramGlobals> GeometryInstance38Type;

typedef Dmrg::MatrixVectorKron<
 Dmrg::ModelBase<
  Dmrg::ModelHelperLocal<
   Dmrg::LeftRightSuper<Dmrg::BasisWithOperators<Dmrg::Operators<Dmrg::Basis<SparseMatrixInstance38Type,CvectorSizeType> >  >,Dmrg::Basis<SparseMatrixInstance38Type,CvectorSizeType> >
  >,
  ParametersDmrgSolverType,
  InputNgType::Readable,
  GeometryInstance38Type
 >
> MatrixVector38Type;

typedef PsimagLite::ChebyshevSolver<PsimagLite::ParametersForSolver<GeometryInstance38Type::RealType>,
	MatrixVector38Type, MatrixVector38Type::VectorType> LanczosSolver38Type;

template void mainLoop4<LanczosSolver38Type,Dmrg::VectorWithOffset<std::complex<RealType> > >
(LanczosSolver38Type::LanczosMatrixType::ModelType::GeometryType&,
const ParametersDmrgSolverType&,
InputNgType::Readable&,
const OperatorOptions&,
PsimagLite::String);


typedef PsimagLite::CrsMatrix<std::complex<RealType> > SparseMatrixInstance39Type;
typedef PsimagLite::Geometry<std::complex<RealType> ,PsimagLite::InputNg<Dmrg::InputCheck>::Readable,Dmrg::ProgramGlobals> GeometryInstance39Type;

typedef Dmrg::MatrixVectorOnTheFly<
 Dmrg::ModelBase<
  Dmrg::ModelHelperLocal<
   Dmrg::LeftRightSuper<Dmrg::BasisWithOperators<Dmrg::Operators<Dmrg::Basis<SparseMatrixInstance39Type,CvectorSizeType> >  >,Dmrg::Basis<SparseMatrixInstance39Type,CvectorSizeType> >
  >,
  ParametersDmrgSolverType,
  InputNgType::Readable,
  GeometryInstance39Type
 >
> MatrixVector39Type;

typedef PsimagLite::ChebyshevSolver<PsimagLite::ParametersForSolver<GeometryInstance39Type::RealType>,
	MatrixVector39Type, MatrixVector39Type::VectorType> LanczosSolver39Type;

template void mainLoop4<LanczosSolver39Type,Dmrg::VectorWithOffsets<std::complex<RealType> > >
(LanczosSolver39Type::LanczosMatrixType::ModelType::GeometryType&,
const ParametersDmrgSolverType&,
InputNgType::Readable&,
const OperatorOptions&,
PsimagLite::String);

