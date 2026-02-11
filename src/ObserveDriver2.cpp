#include "ObserveDriver1.h"

namespace Dmrg {

using QnType                = Qn;
using VectorWithOffset2Type = VectorWithOffset<ComplexType, QnType>;
using VectorWithOffset4Type = VectorWithOffsets<ComplexType, QnType>;

using Geometry1Type = Dmrg::SuperGeometry<RealType, InputNgType::Readable, ProgramGlobals>;

using Basis2Type              = Basis<MySparseMatrixComplex>;
using BasisWithOperators2Type = BasisWithOperators<Basis2Type>;
using LeftRightSuper2Type     = LeftRightSuper<BasisWithOperators2Type, Basis2Type>;

using ModelHelper2Type = ModelHelperLocal<LeftRightSuper2Type>;

using ModelBase5Type
    = ModelBase<ModelHelper2Type, ParametersDmrgSolverType, InputNgType::Readable, Geometry1Type>;

template bool
observeOneFullSweep<VectorWithOffset2Type, ModelBase5Type>(IoInputType&              io,
                                                           const ModelBase5Type&     model,
                                                           const PsimagLite::String& list,
                                                           SizeType                  orbitals);

template bool
observeOneFullSweep<VectorWithOffset4Type, ModelBase5Type>(IoInputType&              io,
                                                           const ModelBase5Type&     model,
                                                           const PsimagLite::String& list,
                                                           SizeType                  orbitals);

}
