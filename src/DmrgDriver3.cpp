// Created automatically by configure.pl
// DO NOT EDIT because file will be overwritten each
// time you run configure.pl with the second argument set to 1
// This file should be commited
#include "DmrgDriver1.h"


typedef PsimagLite::CrsMatrix<RealType> SparseMatrixInstance6Type;
typedef PsimagLite::Geometry<RealType,PsimagLite::InputNg<Dmrg::InputCheck>::Readable,Dmrg::ProgramGlobals> GeometryInstance6Type;

typedef Dmrg::MatrixVectorOnTheFly<
 Dmrg::ModelBase<
  Dmrg::ModelHelperSu2<
   Dmrg::LeftRightSuper<Dmrg::BasisWithOperators<Dmrg::Operators<Dmrg::Basis<SparseMatrixInstance6Type,CvectorSizeType> >  >,Dmrg::Basis<SparseMatrixInstance6Type,CvectorSizeType> >
  >,
  ParametersDmrgSolverType,
  InputNgType::Readable,
  GeometryInstance6Type
 >
> MatrixVector6Type;

typedef PsimagLite::LanczosSolver<PsimagLite::ParametersForSolver<GeometryInstance6Type::RealType>,
	MatrixVector6Type, MatrixVector6Type::VectorType> LanczosSolver6Type;

template void mainLoop4<LanczosSolver6Type,Dmrg::VectorWithOffset<RealType> >
(LanczosSolver6Type::LanczosMatrixType::ModelType::GeometryType&,
const ParametersDmrgSolverType&,
InputNgType::Readable&,
const OperatorOptions&,
PsimagLite::String);

namespace Dmrg {
template<>
bool LinkProductHeisenbergAncillaC<
Dmrg::ModelHelperSu2<
   Dmrg::LeftRightSuper<Dmrg::BasisWithOperators<Dmrg::Operators<Dmrg::Basis<SparseMatrixInstance6Type,CvectorSizeType> >  >,Dmrg::Basis<SparseMatrixInstance6Type,CvectorSizeType> >
  >
>::hot_ = false;

template<>
bool LinkProductHubbardAncillaExtended<
Dmrg::ModelHelperSu2<
   Dmrg::LeftRightSuper<Dmrg::BasisWithOperators<Dmrg::Operators<Dmrg::Basis<SparseMatrixInstance6Type,CvectorSizeType> >  >,Dmrg::Basis<SparseMatrixInstance6Type,CvectorSizeType> >
  >
>::hot_ = false;

template<>
bool LinkProductTjAncillaC2<
Dmrg::ModelHelperSu2<
   Dmrg::LeftRightSuper<Dmrg::BasisWithOperators<Dmrg::Operators<Dmrg::Basis<SparseMatrixInstance6Type,CvectorSizeType> >  >,Dmrg::Basis<SparseMatrixInstance6Type,CvectorSizeType> >
  >
>::hot_ = false;

template<>
bool LinkProdExtendedSuperHubbard1Orb<
Dmrg::ModelHelperSu2<
   Dmrg::LeftRightSuper<Dmrg::BasisWithOperators<Dmrg::Operators<Dmrg::Basis<SparseMatrixInstance6Type,CvectorSizeType> >  >,Dmrg::Basis<SparseMatrixInstance6Type,CvectorSizeType> >
  >
>::hasSpinOrbit_ = false;

template<>
SizeType LinkProductHeisenberg<
Dmrg::ModelHelperSu2<
   Dmrg::LeftRightSuper<Dmrg::BasisWithOperators<Dmrg::Operators<Dmrg::Basis<SparseMatrixInstance6Type,CvectorSizeType> >  >,Dmrg::Basis<SparseMatrixInstance6Type,CvectorSizeType> >
  >
>::terms_ = 2;

} // namespace Dmrg

typedef PsimagLite::CrsMatrix<RealType> SparseMatrixInstance7Type;
typedef PsimagLite::Geometry<RealType,PsimagLite::InputNg<Dmrg::InputCheck>::Readable,Dmrg::ProgramGlobals> GeometryInstance7Type;

typedef Dmrg::MatrixVectorStored<
 Dmrg::ModelBase<
  Dmrg::ModelHelperSu2<
   Dmrg::LeftRightSuper<Dmrg::BasisWithOperators<Dmrg::Operators<Dmrg::Basis<SparseMatrixInstance7Type,CvectorSizeType> >  >,Dmrg::Basis<SparseMatrixInstance7Type,CvectorSizeType> >
  >,
  ParametersDmrgSolverType,
  InputNgType::Readable,
  GeometryInstance7Type
 >
> MatrixVector7Type;

typedef PsimagLite::LanczosSolver<PsimagLite::ParametersForSolver<GeometryInstance7Type::RealType>,
	MatrixVector7Type, MatrixVector7Type::VectorType> LanczosSolver7Type;

template void mainLoop4<LanczosSolver7Type,Dmrg::VectorWithOffset<RealType> >
(LanczosSolver7Type::LanczosMatrixType::ModelType::GeometryType&,
const ParametersDmrgSolverType&,
InputNgType::Readable&,
const OperatorOptions&,
PsimagLite::String);

