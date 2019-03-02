// Created automatically by ./configure.pl
// DO NOT EDIT because file will be overwritten each
// time you run ./configure.pl with the second argument set to 1
// This file should be commited
#include "DmrgDriver1.h"


typedef PsimagLite::CrsMatrix<std::complex<RealType> > SparseMatrixInstance46Type;
typedef PsimagLite::Geometry<std::complex<RealType> ,PsimagLite::InputNg<Dmrg::InputCheck>::Readable,Dmrg::ProgramGlobals> GeometryInstance46Type;

typedef Dmrg::MatrixVectorStored<
 Dmrg::ModelBase<
  Dmrg::ModelHelperSu2<
   Dmrg::LeftRightSuper<Dmrg::BasisWithOperators<Dmrg::Operators<Dmrg::Basis<SparseMatrixInstance46Type> >  >,Dmrg::Basis<SparseMatrixInstance46Type> >
  >,
  ParametersDmrgSolverType,
  InputNgType::Readable,
  GeometryInstance46Type
 >
> MatrixVector46Type;

typedef PsimagLite::ChebyshevSolver<PsimagLite::ParametersForSolver<GeometryInstance46Type::RealType>,
	MatrixVector46Type, MatrixVector46Type::VectorType> LanczosSolver46Type;

template void mainLoop4<LanczosSolver46Type,Dmrg::VectorWithOffsets<std::complex<RealType> , Dmrg::Qn> >
(LanczosSolver46Type::MatrixType::ModelType::GeometryType&,
const ParametersDmrgSolverType&,
InputNgType::Readable&,
const OperatorOptions&);


typedef PsimagLite::CrsMatrix<std::complex<RealType> > SparseMatrixInstance47Type;
typedef PsimagLite::Geometry<std::complex<RealType> ,PsimagLite::InputNg<Dmrg::InputCheck>::Readable,Dmrg::ProgramGlobals> GeometryInstance47Type;

typedef Dmrg::MatrixVectorKron<
 Dmrg::ModelBase<
  Dmrg::ModelHelperSu2<
   Dmrg::LeftRightSuper<Dmrg::BasisWithOperators<Dmrg::Operators<Dmrg::Basis<SparseMatrixInstance47Type> >  >,Dmrg::Basis<SparseMatrixInstance47Type> >
  >,
  ParametersDmrgSolverType,
  InputNgType::Readable,
  GeometryInstance47Type
 >
> MatrixVector47Type;

typedef PsimagLite::ChebyshevSolver<PsimagLite::ParametersForSolver<GeometryInstance47Type::RealType>,
	MatrixVector47Type, MatrixVector47Type::VectorType> LanczosSolver47Type;

template void mainLoop4<LanczosSolver47Type,Dmrg::VectorWithOffsets<std::complex<RealType> , Dmrg::Qn> >
(LanczosSolver47Type::MatrixType::ModelType::GeometryType&,
const ParametersDmrgSolverType&,
InputNgType::Readable&,
const OperatorOptions&);

