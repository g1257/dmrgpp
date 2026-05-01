#include "LanczosDriver1.h"

typedef PsimagLite::Geometry<RealType, InputNgType::Readable, LanczosPlusPlus::LanczosGlobals>
    Geometry1Type;

typedef LanczosPlusPlus::ModelSelector<RealType, Geometry1Type, InputNgType::Readable>
                                          ModelSelector1Type;
typedef ModelSelector1Type::ModelBaseType ModelBase1Type;

typedef ModelBase1Type::BasisBaseType BasisBase1Type;

typedef LanczosPlusPlus::ReflectionSymmetry<Geometry1Type, BasisBase1Type> Symmetry1Type;

template void
mainLoop2<ModelBase1Type, Symmetry1Type>(const ModelBase1Type&            model,
                                         InputNgType::Readable&           io,
                                         LanczosPlusPlus::LanczosOptions& lanczosOptions);
