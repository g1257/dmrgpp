#include "ObserveDriver1.h"

namespace Dmrg {

typedef VectorWithOffset<ComplexType> VectorWithOffset2Type;
typedef VectorWithOffsets<ComplexType> VectorWithOffset4Type;

typedef PsimagLite::Geometry<RealType,
InputNgType::Readable,
ProgramGlobals> Geometry1Type;


typedef Basis<MySparseMatrixComplex, CvectorSizeType> Basis2Type;
typedef Operators<Basis2Type> Operators2Type;
typedef BasisWithOperators<Operators2Type> BasisWithOperators2Type;
typedef LeftRightSuper<BasisWithOperators2Type,Basis2Type> LeftRightSuper2Type;

typedef ModelHelperLocal<LeftRightSuper2Type> ModelHelper2Type;
typedef ModelHelperSu2<LeftRightSuper2Type> ModelHelper4Type;


typedef ModelBase<ModelHelper2Type,
ParametersDmrgSolverType,
InputNgType::Readable,
Geometry1Type> ModelBase5Type;

typedef ModelBase<ModelHelper4Type,
ParametersDmrgSolverType,
InputNgType::Readable,
Geometry1Type> ModelBase6Type;

template bool observeOneFullSweep<VectorWithOffset2Type,ModelBase5Type>(IoInputType& io,
const ModelBase5Type& model,
const PsimagLite::String& list,
bool hasTimeEvolution,
SizeType orbitals);

template bool observeOneFullSweep<VectorWithOffset4Type,ModelBase5Type>(IoInputType& io,
const ModelBase5Type& model,
const PsimagLite::String& list,
bool hasTimeEvolution,
SizeType orbitals);

template bool observeOneFullSweep<VectorWithOffset4Type,ModelBase6Type>(IoInputType& io,
const ModelBase6Type& model,
const PsimagLite::String& list,
bool hasTimeEvolution,
SizeType orbitals);
}
