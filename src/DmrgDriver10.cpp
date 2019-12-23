// Created automatically by ./configure.pl
// DO NOT EDIT because file will be overwritten each
// time you run ./configure.pl with the second argument set to 1
// This file should be commited
#include "DmrgDriver1.h"


typedef PsimagLite::CrsMatrix<std::complex<RealType> > SparseMatrixInstance20Type;
typedef Dmrg::SuperGeometry<std::complex<RealType> ,PsimagLite::InputNg<Dmrg::InputCheck>::Readable,Dmrg::ProgramGlobals> GeometryInstance20Type;

typedef Dmrg::MatrixVectorKron<
 Dmrg::ModelBase<
  Dmrg::ModelHelperLocal<
   Dmrg::LeftRightSuper<Dmrg::BasisWithOperators<Dmrg::Operators<Dmrg::Basis<SparseMatrixInstance20Type> >  >,Dmrg::Basis<SparseMatrixInstance20Type> >
  >,
  ParametersDmrgSolverType,
  InputNgType::Readable,
  GeometryInstance20Type
 >
> MatrixVector20Type;

typedef PsimagLite::ChebyshevSolver<PsimagLite::ParametersForSolver<GeometryInstance20Type::RealType>,
	MatrixVector20Type, MatrixVector20Type::VectorType> LanczosSolver20Type;

template void mainLoop4<LanczosSolver20Type,Dmrg::VectorWithOffset<std::complex<RealType> , Dmrg::Qn> >
(LanczosSolver20Type::MatrixType::ModelType::SuperGeometryType&,
const ParametersDmrgSolverType&,
InputNgType::Readable&,
const OperatorOptions&);


typedef PsimagLite::CrsMatrix<std::complex<RealType> > SparseMatrixInstance21Type;
typedef Dmrg::SuperGeometry<std::complex<RealType> ,PsimagLite::InputNg<Dmrg::InputCheck>::Readable,Dmrg::ProgramGlobals> GeometryInstance21Type;

typedef Dmrg::MatrixVectorOnTheFly<
 Dmrg::ModelBase<
  Dmrg::ModelHelperLocal<
   Dmrg::LeftRightSuper<Dmrg::BasisWithOperators<Dmrg::Operators<Dmrg::Basis<SparseMatrixInstance21Type> >  >,Dmrg::Basis<SparseMatrixInstance21Type> >
  >,
  ParametersDmrgSolverType,
  InputNgType::Readable,
  GeometryInstance21Type
 >
> MatrixVector21Type;

typedef PsimagLite::ChebyshevSolver<PsimagLite::ParametersForSolver<GeometryInstance21Type::RealType>,
	MatrixVector21Type, MatrixVector21Type::VectorType> LanczosSolver21Type;

template void mainLoop4<LanczosSolver21Type,Dmrg::VectorWithOffsets<std::complex<RealType> , Dmrg::Qn> >
(LanczosSolver21Type::MatrixType::ModelType::SuperGeometryType&,
const ParametersDmrgSolverType&,
InputNgType::Readable&,
const OperatorOptions&);

