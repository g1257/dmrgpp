// Created automatically by ./newconfigure.pl
// DO NOT EDIT because file will be overwritten each
// time you run ./newconfigure.pl with the second argument set to 1
// This file should be commited
#include "DmrgDriver1.h"


typedef PsimagLite::CrsMatrix<RealType> SparseMatrixInstance14Type;
typedef PsimagLite::Geometry<RealType,PsimagLite::InputNg<Dmrg::InputCheck>::Readable,Dmrg::ProgramGlobals> GeometryInstance14Type;

typedef Dmrg::MatrixVectorKron<
 Dmrg::ModelBase<
  Dmrg::ModelHelperLocal<
   Dmrg::LeftRightSuper<Dmrg::BasisWithOperators<Dmrg::Operators<Dmrg::Basis<SparseMatrixInstance14Type> >  >,Dmrg::Basis<SparseMatrixInstance14Type> >
  >,
  ParametersDmrgSolverType,
  InputNgType::Readable,
  GeometryInstance14Type
 >
> MatrixVector14Type;

typedef PsimagLite::ChebyshevSolver<PsimagLite::ParametersForSolver<GeometryInstance14Type::RealType>,
	MatrixVector14Type, MatrixVector14Type::VectorType> LanczosSolver14Type;

template void mainLoop4<LanczosSolver14Type,Dmrg::VectorWithOffset<RealType, Dmrg::Qn> >
(LanczosSolver14Type::LanczosMatrixType::ModelType::GeometryType&,
const ParametersDmrgSolverType&,
InputNgType::Readable&,
const OperatorOptions&,
PsimagLite::String);


typedef PsimagLite::CrsMatrix<RealType> SparseMatrixInstance15Type;
typedef PsimagLite::Geometry<RealType,PsimagLite::InputNg<Dmrg::InputCheck>::Readable,Dmrg::ProgramGlobals> GeometryInstance15Type;

typedef Dmrg::MatrixVectorOnTheFly<
 Dmrg::ModelBase<
  Dmrg::ModelHelperLocal<
   Dmrg::LeftRightSuper<Dmrg::BasisWithOperators<Dmrg::Operators<Dmrg::Basis<SparseMatrixInstance15Type> >  >,Dmrg::Basis<SparseMatrixInstance15Type> >
  >,
  ParametersDmrgSolverType,
  InputNgType::Readable,
  GeometryInstance15Type
 >
> MatrixVector15Type;

typedef PsimagLite::ChebyshevSolver<PsimagLite::ParametersForSolver<GeometryInstance15Type::RealType>,
	MatrixVector15Type, MatrixVector15Type::VectorType> LanczosSolver15Type;

template void mainLoop4<LanczosSolver15Type,Dmrg::VectorWithOffsets<RealType, Dmrg::Qn> >
(LanczosSolver15Type::LanczosMatrixType::ModelType::GeometryType&,
const ParametersDmrgSolverType&,
InputNgType::Readable&,
const OperatorOptions&,
PsimagLite::String);

