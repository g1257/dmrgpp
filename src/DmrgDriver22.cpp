// Created automatically by ./configure.pl
// DO NOT EDIT because file will be overwritten each
// time you run ./configure.pl with the second argument set to 1
// This file should be commited
#include "DmrgDriver1.h"


typedef PsimagLite::CrsMatrix<std::complex<RealType> > SparseMatrixInstance44Type;
typedef PsimagLite::Geometry<std::complex<RealType> ,PsimagLite::InputNg<Dmrg::InputCheck>::Readable,Dmrg::ProgramGlobals> GeometryInstance44Type;

typedef Dmrg::MatrixVectorKron<
 Dmrg::ModelBase<
  Dmrg::ModelHelperSu2<
   Dmrg::LeftRightSuper<Dmrg::BasisWithOperators<Dmrg::Operators<Dmrg::Basis<SparseMatrixInstance44Type> >  >,Dmrg::Basis<SparseMatrixInstance44Type> >
  >,
  ParametersDmrgSolverType,
  InputNgType::Readable,
  GeometryInstance44Type
 >
> MatrixVector44Type;

typedef PsimagLite::ChebyshevSolver<PsimagLite::ParametersForSolver<GeometryInstance44Type::RealType>,
	MatrixVector44Type, MatrixVector44Type::VectorType> LanczosSolver44Type;

template void mainLoop4<LanczosSolver44Type,Dmrg::VectorWithOffset<std::complex<RealType> , Dmrg::Qn> >
(LanczosSolver44Type::MatrixType::ModelType::GeometryType&,
const ParametersDmrgSolverType&,
InputNgType::Readable&,
const OperatorOptions&,
PsimagLite::String);


typedef PsimagLite::CrsMatrix<std::complex<RealType> > SparseMatrixInstance45Type;
typedef PsimagLite::Geometry<std::complex<RealType> ,PsimagLite::InputNg<Dmrg::InputCheck>::Readable,Dmrg::ProgramGlobals> GeometryInstance45Type;

typedef Dmrg::MatrixVectorOnTheFly<
 Dmrg::ModelBase<
  Dmrg::ModelHelperSu2<
   Dmrg::LeftRightSuper<Dmrg::BasisWithOperators<Dmrg::Operators<Dmrg::Basis<SparseMatrixInstance45Type> >  >,Dmrg::Basis<SparseMatrixInstance45Type> >
  >,
  ParametersDmrgSolverType,
  InputNgType::Readable,
  GeometryInstance45Type
 >
> MatrixVector45Type;

typedef PsimagLite::ChebyshevSolver<PsimagLite::ParametersForSolver<GeometryInstance45Type::RealType>,
	MatrixVector45Type, MatrixVector45Type::VectorType> LanczosSolver45Type;

template void mainLoop4<LanczosSolver45Type,Dmrg::VectorWithOffsets<std::complex<RealType> , Dmrg::Qn> >
(LanczosSolver45Type::MatrixType::ModelType::GeometryType&,
const ParametersDmrgSolverType&,
InputNgType::Readable&,
const OperatorOptions&,
PsimagLite::String);

