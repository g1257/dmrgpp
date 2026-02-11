#include "ObserveDriver1.h"

namespace Dmrg {

using QnType                = Qn;
using VectorWithOffset1Type = VectorWithOffset<RealType, QnType>;
using VectorWithOffset2Type = VectorWithOffset<ComplexType, QnType>;
using VectorWithOffset3Type = VectorWithOffsets<RealType, QnType>;
using VectorWithOffset4Type = VectorWithOffsets<ComplexType, QnType>;

using Geometry1Type = Dmrg::SuperGeometry<RealType, InputNgType::Readable, ProgramGlobals>;

using Geometry2Type = Dmrg::SuperGeometry<ComplexType, InputNgType::Readable, ProgramGlobals>;

using Basis1Type              = Basis<MySparseMatrixReal>;
using BasisWithOperators1Type = BasisWithOperators<Basis1Type>;
using LeftRightSuper1Type     = LeftRightSuper<BasisWithOperators1Type, Basis1Type>;

using Basis2Type              = Basis<MySparseMatrixComplex>;
using BasisWithOperators2Type = BasisWithOperators<Basis2Type>;
using LeftRightSuper2Type     = LeftRightSuper<BasisWithOperators2Type, Basis2Type>;

}
