#include "ObserveDriver1.h"

namespace Dmrg {

typedef Qn QnType;
typedef VectorWithOffset<ComplexType, QnType> VectorWithOffset2Type;
typedef VectorWithOffsets<ComplexType, QnType> VectorWithOffset4Type;

typedef Dmrg::SuperGeometry<RealType,
InputNgType::Readable,
ProgramGlobals> Geometry1Type;


typedef Basis<MySparseMatrixComplex> Basis2Type;
typedef Operators<Basis2Type> Operators2Type;
typedef BasisWithOperators<Operators2Type> BasisWithOperators2Type;
typedef LeftRightSuper<BasisWithOperators2Type,Basis2Type> LeftRightSuper2Type;

typedef ModelHelperLocal<LeftRightSuper2Type> ModelHelper2Type;

typedef ModelBase<ModelHelper2Type,
ParametersDmrgSolverType,
InputNgType::Readable,
Geometry1Type> ModelBase5Type;

template bool observeOneFullSweep<VectorWithOffset2Type,ModelBase5Type>(IoInputType& io,
const ModelBase5Type& model,
const PsimagLite::String& list,
SizeType orbitals);

template bool observeOneFullSweep<VectorWithOffset4Type,ModelBase5Type>(IoInputType& io,
const ModelBase5Type& model,
const PsimagLite::String& list,
SizeType orbitals);

}
