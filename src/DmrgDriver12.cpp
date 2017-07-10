// Created automatically by configure.pl
// DO NOT EDIT because file will be overwritten each
// time you run configure.pl with the second argument set to 1
// This file should be commited
#include "DmrgDriver1.h"


typedef PsimagLite::CrsMatrix<std::complex<RealType> > SparseMatrixInstance24Type;
typedef PsimagLite::Geometry<std::complex<RealType> ,PsimagLite::InputNg<Dmrg::InputCheck>::Readable,Dmrg::ProgramGlobals> GeometryInstance24Type;

typedef Dmrg::MatrixVectorOnTheFly<
 Dmrg::ModelBase<
  Dmrg::ModelHelperLocal<
   Dmrg::LeftRightSuper<Dmrg::BasisWithOperators<Dmrg::Operators<Dmrg::Basis<SparseMatrixInstance24Type,CvectorSizeType> >  >,Dmrg::Basis<SparseMatrixInstance24Type,CvectorSizeType> >
  >,
  ParametersDmrgSolverType,
  InputNgType::Readable,
  GeometryInstance24Type
 >
> MatrixVector24Type;

typedef PsimagLite::LanczosSolver<PsimagLite::ParametersForSolver<GeometryInstance24Type::RealType>,
	MatrixVector24Type, MatrixVector24Type::VectorType> LanczosSolver24Type;

template void mainLoop4<LanczosSolver24Type,Dmrg::VectorWithOffset<std::complex<RealType> > >
(LanczosSolver24Type::LanczosMatrixType::ModelType::GeometryType&,
const ParametersDmrgSolverType&,
InputNgType::Readable&,
const OperatorOptions&,
PsimagLite::String);

namespace Dmrg {
template<>
bool LinkProductHeisenbergAncillaC<
Dmrg::ModelHelperLocal<
   Dmrg::LeftRightSuper<Dmrg::BasisWithOperators<Dmrg::Operators<Dmrg::Basis<SparseMatrixInstance24Type,CvectorSizeType> >  >,Dmrg::Basis<SparseMatrixInstance24Type,CvectorSizeType> >
  >
>::hot_ = false;

template<>
bool LinkProductHubbardAncillaExtended<
Dmrg::ModelHelperLocal<
   Dmrg::LeftRightSuper<Dmrg::BasisWithOperators<Dmrg::Operators<Dmrg::Basis<SparseMatrixInstance24Type,CvectorSizeType> >  >,Dmrg::Basis<SparseMatrixInstance24Type,CvectorSizeType> >
  >
>::hot_ = false;

template<>
bool LinkProductTjAncillaC2<
Dmrg::ModelHelperLocal<
   Dmrg::LeftRightSuper<Dmrg::BasisWithOperators<Dmrg::Operators<Dmrg::Basis<SparseMatrixInstance24Type,CvectorSizeType> >  >,Dmrg::Basis<SparseMatrixInstance24Type,CvectorSizeType> >
  >
>::hot_ = false;

template<>
bool LinkProdExtendedSuperHubbard1Orb<
Dmrg::ModelHelperLocal<
   Dmrg::LeftRightSuper<Dmrg::BasisWithOperators<Dmrg::Operators<Dmrg::Basis<SparseMatrixInstance24Type,CvectorSizeType> >  >,Dmrg::Basis<SparseMatrixInstance24Type,CvectorSizeType> >
  >
>::hasSpinOrbit_ = false;

template<>
SizeType LinkProductHeisenberg<
Dmrg::ModelHelperLocal<
   Dmrg::LeftRightSuper<Dmrg::BasisWithOperators<Dmrg::Operators<Dmrg::Basis<SparseMatrixInstance24Type,CvectorSizeType> >  >,Dmrg::Basis<SparseMatrixInstance24Type,CvectorSizeType> >
  >
>::terms_ = 2;

} // namespace Dmrg

typedef PsimagLite::CrsMatrix<std::complex<RealType> > SparseMatrixInstance25Type;
typedef PsimagLite::Geometry<std::complex<RealType> ,PsimagLite::InputNg<Dmrg::InputCheck>::Readable,Dmrg::ProgramGlobals> GeometryInstance25Type;

typedef Dmrg::MatrixVectorStored<
 Dmrg::ModelBase<
  Dmrg::ModelHelperLocal<
   Dmrg::LeftRightSuper<Dmrg::BasisWithOperators<Dmrg::Operators<Dmrg::Basis<SparseMatrixInstance25Type,CvectorSizeType> >  >,Dmrg::Basis<SparseMatrixInstance25Type,CvectorSizeType> >
  >,
  ParametersDmrgSolverType,
  InputNgType::Readable,
  GeometryInstance25Type
 >
> MatrixVector25Type;

typedef PsimagLite::LanczosSolver<PsimagLite::ParametersForSolver<GeometryInstance25Type::RealType>,
	MatrixVector25Type, MatrixVector25Type::VectorType> LanczosSolver25Type;

template void mainLoop4<LanczosSolver25Type,Dmrg::VectorWithOffset<std::complex<RealType> > >
(LanczosSolver25Type::LanczosMatrixType::ModelType::GeometryType&,
const ParametersDmrgSolverType&,
InputNgType::Readable&,
const OperatorOptions&,
PsimagLite::String);

