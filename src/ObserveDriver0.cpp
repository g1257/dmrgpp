#include "ObserveDriver1.h"

namespace Dmrg {

typedef Qn QnType;
typedef VectorWithOffset<RealType, QnType> VectorWithOffset1Type;
typedef VectorWithOffset<ComplexType, QnType> VectorWithOffset2Type;
typedef VectorWithOffsets<RealType, QnType> VectorWithOffset3Type;
typedef VectorWithOffsets<ComplexType, QnType> VectorWithOffset4Type;

typedef Dmrg::SuperGeometry<RealType,
InputNgType::Readable,
ProgramGlobals> Geometry1Type;

typedef Dmrg::SuperGeometry<ComplexType,
InputNgType::Readable,
ProgramGlobals> Geometry2Type;

typedef Basis<MySparseMatrixReal> Basis1Type;
typedef BasisWithOperators<Basis1Type> BasisWithOperators1Type;
typedef LeftRightSuper<BasisWithOperators1Type,Basis1Type> LeftRightSuper1Type;

typedef Basis<MySparseMatrixComplex> Basis2Type;
typedef BasisWithOperators<Basis2Type> BasisWithOperators2Type;
typedef LeftRightSuper<BasisWithOperators2Type,Basis2Type> LeftRightSuper2Type;

typedef ModelHelperLocal<LeftRightSuper1Type> ModelHelper1Type;
typedef ModelHelperLocal<LeftRightSuper2Type> ModelHelper2Type;

typedef ModelBase<ModelHelper1Type,
ParametersDmrgSolverType,
InputNgType::Readable,
Geometry1Type> ModelBase1Type;

typedef ModelBase<ModelHelper2Type,
ParametersDmrgSolverType,
InputNgType::Readable,
Geometry2Type> ModelBase2Type;


template bool observeOneFullSweep<VectorWithOffset1Type,ModelBase1Type>(IoInputType& io,
const ModelBase1Type& model,
const PsimagLite::String& list,
SizeType orbitals);

template bool observeOneFullSweep<VectorWithOffset2Type,ModelBase2Type>(IoInputType& io,
const ModelBase2Type& model,
const PsimagLite::String& list,
SizeType orbitals);

template bool observeOneFullSweep<VectorWithOffset3Type,ModelBase1Type>(IoInputType& io,
const ModelBase1Type& model,
const PsimagLite::String& list,
SizeType orbitals);

template bool observeOneFullSweep<VectorWithOffset4Type,ModelBase2Type>(IoInputType& io,
const ModelBase2Type& model,
const PsimagLite::String& list,
SizeType orbitals);
}
