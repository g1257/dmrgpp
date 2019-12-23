// Created automatically by ./configure.pl
// DO NOT EDIT because file will be overwritten each
// time you run ./configure.pl with the second argument set to 1
// This file should be commited
#include "DmrgDriver1.h"


typedef PsimagLite::CrsMatrix<std::complex<RealType> > SparseMatrixInstance22Type;
typedef Dmrg::SuperGeometry<std::complex<RealType> ,PsimagLite::InputNg<Dmrg::InputCheck>::Readable,Dmrg::ProgramGlobals> GeometryInstance22Type;

typedef Dmrg::MatrixVectorStored<
 Dmrg::ModelBase<
  Dmrg::ModelHelperLocal<
   Dmrg::LeftRightSuper<Dmrg::BasisWithOperators<Dmrg::Operators<Dmrg::Basis<SparseMatrixInstance22Type> >  >,Dmrg::Basis<SparseMatrixInstance22Type> >
  >,
  ParametersDmrgSolverType,
  InputNgType::Readable,
  GeometryInstance22Type
 >
> MatrixVector22Type;

typedef PsimagLite::ChebyshevSolver<PsimagLite::ParametersForSolver<GeometryInstance22Type::RealType>,
	MatrixVector22Type, MatrixVector22Type::VectorType> LanczosSolver22Type;

template void mainLoop4<LanczosSolver22Type,Dmrg::VectorWithOffsets<std::complex<RealType> , Dmrg::Qn> >
(LanczosSolver22Type::MatrixType::ModelType::SuperGeometryType&,
const ParametersDmrgSolverType&,
InputNgType::Readable&,
const OperatorOptions&);


typedef PsimagLite::CrsMatrix<std::complex<RealType> > SparseMatrixInstance23Type;
typedef Dmrg::SuperGeometry<std::complex<RealType> ,PsimagLite::InputNg<Dmrg::InputCheck>::Readable,Dmrg::ProgramGlobals> GeometryInstance23Type;

typedef Dmrg::MatrixVectorKron<
 Dmrg::ModelBase<
  Dmrg::ModelHelperLocal<
   Dmrg::LeftRightSuper<Dmrg::BasisWithOperators<Dmrg::Operators<Dmrg::Basis<SparseMatrixInstance23Type> >  >,Dmrg::Basis<SparseMatrixInstance23Type> >
  >,
  ParametersDmrgSolverType,
  InputNgType::Readable,
  GeometryInstance23Type
 >
> MatrixVector23Type;

typedef PsimagLite::ChebyshevSolver<PsimagLite::ParametersForSolver<GeometryInstance23Type::RealType>,
	MatrixVector23Type, MatrixVector23Type::VectorType> LanczosSolver23Type;

template void mainLoop4<LanczosSolver23Type,Dmrg::VectorWithOffsets<std::complex<RealType> , Dmrg::Qn> >
(LanczosSolver23Type::MatrixType::ModelType::SuperGeometryType&,
const ParametersDmrgSolverType&,
InputNgType::Readable&,
const OperatorOptions&);

