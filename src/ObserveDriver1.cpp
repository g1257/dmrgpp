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
typedef BasisWithOperators<Basis1Type, 1> BasisWithOperators1Type;
typedef LeftRightSuper<BasisWithOperators1Type,Basis1Type> LeftRightSuper1Type;

typedef Basis<MySparseMatrixComplex> Basis2Type;
typedef BasisWithOperators<Basis2Type, 1> BasisWithOperators2Type;
typedef LeftRightSuper<BasisWithOperators2Type,Basis2Type> LeftRightSuper2Type;

}
