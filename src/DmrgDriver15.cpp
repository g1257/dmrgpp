// Created automatically by configure.pl
// DO NOT EDIT because file will be overwritten each
// time you run configure.pl with the second argument set to 1
// This file should be commited
#include "DmrgDriver1.h"


typedef PsimagLite::CrsMatrix<std::complex<RealType> > SparseMatrixInstance30Type;
typedef PsimagLite::Geometry<std::complex<RealType> ,PsimagLite::InputNg<Dmrg::InputCheck>::Readable,Dmrg::ProgramGlobals> GeometryInstance30Type;

typedef Dmrg::MatrixVectorOnTheFly<
 Dmrg::ModelBase<
  Dmrg::ModelHelperSu2<
   Dmrg::LeftRightSuper<Dmrg::BasisWithOperators<Dmrg::Operators<Dmrg::Basis<SparseMatrixInstance30Type,CvectorSizeType> >  >,Dmrg::Basis<SparseMatrixInstance30Type,CvectorSizeType> >
  >,
  ParametersDmrgSolverType,
  InputNgType::Readable,
  GeometryInstance30Type
 >
> MatrixVector30Type;

typedef PsimagLite::LanczosSolver<PsimagLite::ParametersForSolver<GeometryInstance30Type::RealType>,
	MatrixVector30Type, MatrixVector30Type::VectorType> LanczosSolver30Type;

template void mainLoop4<LanczosSolver30Type,Dmrg::VectorWithOffset<std::complex<RealType> > >
(LanczosSolver30Type::LanczosMatrixType::ModelType::GeometryType&,
const ParametersDmrgSolverType&,
InputNgType::Readable&,
const OperatorOptions&,
PsimagLite::String);

namespace Dmrg {
template<>
bool LinkProductHeisenbergAncillaC<
Dmrg::ModelHelperSu2<
   Dmrg::LeftRightSuper<Dmrg::BasisWithOperators<Dmrg::Operators<Dmrg::Basis<SparseMatrixInstance30Type,CvectorSizeType> >  >,Dmrg::Basis<SparseMatrixInstance30Type,CvectorSizeType> >
  >
>::hot_ = false;

template<>
bool LinkProductHubbardAncillaExtended<
Dmrg::ModelHelperSu2<
   Dmrg::LeftRightSuper<Dmrg::BasisWithOperators<Dmrg::Operators<Dmrg::Basis<SparseMatrixInstance30Type,CvectorSizeType> >  >,Dmrg::Basis<SparseMatrixInstance30Type,CvectorSizeType> >
  >
>::hot_ = false;

template<>
bool LinkProductTjAncillaC2<
Dmrg::ModelHelperSu2<
   Dmrg::LeftRightSuper<Dmrg::BasisWithOperators<Dmrg::Operators<Dmrg::Basis<SparseMatrixInstance30Type,CvectorSizeType> >  >,Dmrg::Basis<SparseMatrixInstance30Type,CvectorSizeType> >
  >
>::hot_ = false;

template<>
bool LinkProdExtendedSuperHubbard1Orb<
Dmrg::ModelHelperSu2<
   Dmrg::LeftRightSuper<Dmrg::BasisWithOperators<Dmrg::Operators<Dmrg::Basis<SparseMatrixInstance30Type,CvectorSizeType> >  >,Dmrg::Basis<SparseMatrixInstance30Type,CvectorSizeType> >
  >
>::hasSpinOrbit_ = false;
}

typedef PsimagLite::CrsMatrix<std::complex<RealType> > SparseMatrixInstance31Type;
typedef PsimagLite::Geometry<std::complex<RealType> ,PsimagLite::InputNg<Dmrg::InputCheck>::Readable,Dmrg::ProgramGlobals> GeometryInstance31Type;

typedef Dmrg::MatrixVectorStored<
 Dmrg::ModelBase<
  Dmrg::ModelHelperSu2<
   Dmrg::LeftRightSuper<Dmrg::BasisWithOperators<Dmrg::Operators<Dmrg::Basis<SparseMatrixInstance31Type,CvectorSizeType> >  >,Dmrg::Basis<SparseMatrixInstance31Type,CvectorSizeType> >
  >,
  ParametersDmrgSolverType,
  InputNgType::Readable,
  GeometryInstance31Type
 >
> MatrixVector31Type;

typedef PsimagLite::LanczosSolver<PsimagLite::ParametersForSolver<GeometryInstance31Type::RealType>,
	MatrixVector31Type, MatrixVector31Type::VectorType> LanczosSolver31Type;

template void mainLoop4<LanczosSolver31Type,Dmrg::VectorWithOffset<std::complex<RealType> > >
(LanczosSolver31Type::LanczosMatrixType::ModelType::GeometryType&,
const ParametersDmrgSolverType&,
InputNgType::Readable&,
const OperatorOptions&,
PsimagLite::String);

