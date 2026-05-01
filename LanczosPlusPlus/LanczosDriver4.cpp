#include "LanczosDriver1.h"

typedef PsimagLite::
    Geometry<std::complex<RealType>, InputNgType::Readable, LanczosPlusPlus::LanczosGlobals>
        Geometry4Type;

typedef LanczosPlusPlus::ModelSelector<std::complex<RealType>, Geometry4Type, InputNgType::Readable>
                                          ModelSelector4Type;
typedef ModelSelector4Type::ModelBaseType ModelBase4Type;

typedef ModelBase4Type::BasisBaseType BasisBase4Type;

typedef LanczosPlusPlus::ReflectionSymmetry<Geometry4Type, BasisBase4Type> Symmetry4Type;

template void
mainLoop2<ModelBase4Type, Symmetry4Type>(const ModelBase4Type&            model,
                                         InputNgType::Readable&           io,
                                         LanczosPlusPlus::LanczosOptions& lanczosOptions);
