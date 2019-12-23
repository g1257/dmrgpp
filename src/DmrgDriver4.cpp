// Created automatically by ./configure.pl
// DO NOT EDIT because file will be overwritten each
// time you run ./configure.pl with the second argument set to 1
// This file should be commited
#include "DmrgDriver1.h"


typedef PsimagLite::CrsMatrix<RealType> SparseMatrixInstance8Type;
typedef Dmrg::SuperGeometry<RealType,PsimagLite::InputNg<Dmrg::InputCheck>::Readable,Dmrg::ProgramGlobals> GeometryInstance8Type;

typedef Dmrg::MatrixVectorKron<
 Dmrg::ModelBase<
  Dmrg::ModelHelperLocal<
   Dmrg::LeftRightSuper<Dmrg::BasisWithOperators<Dmrg::Operators<Dmrg::Basis<SparseMatrixInstance8Type> >  >,Dmrg::Basis<SparseMatrixInstance8Type> >
  >,
  ParametersDmrgSolverType,
  InputNgType::Readable,
  GeometryInstance8Type
 >
> MatrixVector8Type;

typedef PsimagLite::ChebyshevSolver<PsimagLite::ParametersForSolver<GeometryInstance8Type::RealType>,
	MatrixVector8Type, MatrixVector8Type::VectorType> LanczosSolver8Type;

template void mainLoop4<LanczosSolver8Type,Dmrg::VectorWithOffset<RealType, Dmrg::Qn> >
(LanczosSolver8Type::MatrixType::ModelType::SuperGeometryType&,
const ParametersDmrgSolverType&,
InputNgType::Readable&,
const OperatorOptions&);


typedef PsimagLite::CrsMatrix<RealType> SparseMatrixInstance9Type;
typedef Dmrg::SuperGeometry<RealType,PsimagLite::InputNg<Dmrg::InputCheck>::Readable,Dmrg::ProgramGlobals> GeometryInstance9Type;

typedef Dmrg::MatrixVectorOnTheFly<
 Dmrg::ModelBase<
  Dmrg::ModelHelperLocal<
   Dmrg::LeftRightSuper<Dmrg::BasisWithOperators<Dmrg::Operators<Dmrg::Basis<SparseMatrixInstance9Type> >  >,Dmrg::Basis<SparseMatrixInstance9Type> >
  >,
  ParametersDmrgSolverType,
  InputNgType::Readable,
  GeometryInstance9Type
 >
> MatrixVector9Type;

typedef PsimagLite::ChebyshevSolver<PsimagLite::ParametersForSolver<GeometryInstance9Type::RealType>,
	MatrixVector9Type, MatrixVector9Type::VectorType> LanczosSolver9Type;

template void mainLoop4<LanczosSolver9Type,Dmrg::VectorWithOffsets<RealType, Dmrg::Qn> >
(LanczosSolver9Type::MatrixType::ModelType::SuperGeometryType&,
const ParametersDmrgSolverType&,
InputNgType::Readable&,
const OperatorOptions&);

