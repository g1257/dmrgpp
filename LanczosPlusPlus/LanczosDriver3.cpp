#include "LanczosDriver1.h"

typedef PsimagLite::
    Geometry<std::complex<RealType>, InputNgType::Readable, LanczosPlusPlus::LanczosGlobals>
        Geometry3Type;

typedef LanczosPlusPlus::ModelSelector<std::complex<RealType>, Geometry3Type, InputNgType::Readable>
                                          ModelSelector3Type;
typedef ModelSelector3Type::ModelBaseType ModelBase3Type;

typedef ModelBase3Type::BasisBaseType BasisBase3Type;

typedef LanczosPlusPlus::TranslationSymmetry<Geometry3Type, BasisBase3Type> Symmetry3Type;

template void
mainLoop2<ModelBase3Type, Symmetry3Type>(const ModelBase3Type&            model,
                                         InputNgType::Readable&           io,
                                         LanczosPlusPlus::LanczosOptions& lanczosOptions);
