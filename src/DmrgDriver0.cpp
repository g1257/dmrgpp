// Created automatically by configure.pl
// DO NOT EDIT because file will be overwritten each
// time you run configure.pl with the second argument set to 1
// This file should be commited
#include "DmrgDriver1.h"


typedef PsimagLite::CrsMatrix<RealType> SparseMatrixInstance0Type;
typedef PsimagLite::Geometry<RealType,PsimagLite::InputNg<Dmrg::InputCheck>::Readable,Dmrg::ProgramGlobals> GeometryInstance0Type;

typedef Dmrg::MatrixVectorOnTheFly<
 Dmrg::ModelBase<
  Dmrg::ModelHelperLocal<
   Dmrg::LeftRightSuper<Dmrg::BasisWithOperators<Dmrg::Operators<Dmrg::Basis<SparseMatrixInstance0Type,CvectorSizeType> >  >,Dmrg::Basis<SparseMatrixInstance0Type,CvectorSizeType> >
  >,
  ParametersDmrgSolverType,
  InputNgType::Readable,
  GeometryInstance0Type
 >
> MatrixVector0Type;

typedef PsimagLite::LanczosSolver<PsimagLite::ParametersForSolver<GeometryInstance0Type::RealType>,
	MatrixVector0Type, MatrixVector0Type::VectorType> LanczosSolver0Type;

template void mainLoop4<LanczosSolver0Type,Dmrg::VectorWithOffset<RealType> >
(LanczosSolver0Type::LanczosMatrixType::ModelType::GeometryType&,
const ParametersDmrgSolverType&,
InputNgType::Readable&,
const OperatorOptions&,
PsimagLite::String);


typedef PsimagLite::CrsMatrix<RealType> SparseMatrixInstance1Type;
typedef PsimagLite::Geometry<RealType,PsimagLite::InputNg<Dmrg::InputCheck>::Readable,Dmrg::ProgramGlobals> GeometryInstance1Type;

typedef Dmrg::MatrixVectorStored<
 Dmrg::ModelBase<
  Dmrg::ModelHelperLocal<
   Dmrg::LeftRightSuper<Dmrg::BasisWithOperators<Dmrg::Operators<Dmrg::Basis<SparseMatrixInstance1Type,CvectorSizeType> >  >,Dmrg::Basis<SparseMatrixInstance1Type,CvectorSizeType> >
  >,
  ParametersDmrgSolverType,
  InputNgType::Readable,
  GeometryInstance1Type
 >
> MatrixVector1Type;

typedef PsimagLite::LanczosSolver<PsimagLite::ParametersForSolver<GeometryInstance1Type::RealType>,
	MatrixVector1Type, MatrixVector1Type::VectorType> LanczosSolver1Type;

template void mainLoop4<LanczosSolver1Type,Dmrg::VectorWithOffset<RealType> >
(LanczosSolver1Type::LanczosMatrixType::ModelType::GeometryType&,
const ParametersDmrgSolverType&,
InputNgType::Readable&,
const OperatorOptions&,
PsimagLite::String);

