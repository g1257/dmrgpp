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

using ModelHelper1Type = ModelHelperLocal<LeftRightSuper1Type>;
using ModelHelper2Type = ModelHelperLocal<LeftRightSuper2Type>;

using ModelBase1Type
    = ModelBase<ModelHelper1Type, ParametersDmrgSolverType, InputNgType::Readable, Geometry1Type>;

using ModelBase2Type
    = ModelBase<ModelHelper2Type, ParametersDmrgSolverType, InputNgType::Readable, Geometry2Type>;

template bool
observeOneFullSweep<VectorWithOffset1Type, ModelBase1Type>(IoInputType&              io,
                                                           const ModelBase1Type&     model,
                                                           const PsimagLite::String& list,
                                                           SizeType                  orbitals);

template bool
observeOneFullSweep<VectorWithOffset2Type, ModelBase2Type>(IoInputType&              io,
                                                           const ModelBase2Type&     model,
                                                           const PsimagLite::String& list,
                                                           SizeType                  orbitals);

template bool
observeOneFullSweep<VectorWithOffset3Type, ModelBase1Type>(IoInputType&              io,
                                                           const ModelBase1Type&     model,
                                                           const PsimagLite::String& list,
                                                           SizeType                  orbitals);

template bool
observeOneFullSweep<VectorWithOffset4Type, ModelBase2Type>(IoInputType&              io,
                                                           const ModelBase2Type&     model,
                                                           const PsimagLite::String& list,
                                                           SizeType                  orbitals);
}
