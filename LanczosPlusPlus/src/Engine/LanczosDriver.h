#ifndef LANCZOSDRIVER_H
#define LANCZOSDRIVER_H
#include "../Version.h"
#include "DefaultSymmetry.h"
#include "Engine.h"
#include "InputCheck.h"
#include "InternalProductOnTheFly.h"
#include "InternalProductStored.h"
#include "LanczosGlobals.h"
#include "LanczosOptions.h"
#include "ModelSelector.h"
#include "ReducedDensityMatrix.h"
#include "ReflectionSymmetry.h"
#include "TranslationSymmetry.h"
#include <PsimagLite/AllocatorCpu.h>
#include <PsimagLite/Concurrency.h>
#include <PsimagLite/ContinuedFraction.h> // in PsimagLite
#include <PsimagLite/ContinuedFractionCollection.h> // in PsimagLite
#include <PsimagLite/Geometry/Geometry.h>
#include <PsimagLite/InputNg.h> // in PsimagLite
#include <PsimagLite/PsimagLite.h>
#include <PsimagLite/Version.h>
#include <cstdlib>
#include <getopt.h>
#include <unistd.h>

#ifndef USE_FLOAT
typedef double RealType;
#else
typedef float RealType;
#endif

typedef PsimagLite::Concurrency                          ConcurrencyType;
typedef PsimagLite::InputNg<LanczosPlusPlus::InputCheck> InputNgType;

template <typename ModelType> SizeType maxOrbitals(const ModelType& model);

template <typename EngineType>
void extendedStatic(PsimagLite::String                   manypoint,
                    const EngineType&                    engine,
                    const typename EngineType::PairType& braAndKet);

template <typename ModelType,
          typename SpecialSymmetryType,
          template <typename, typename> class InternalProductTemplate>
void mainLoop3(const ModelType&                 model,
               InputNgType::Readable&           io,
               LanczosPlusPlus::LanczosOptions& lanczosOptions);

template <typename ModelType, typename SpecialSymmetryType>
void mainLoop2(const ModelType&                 model,
               InputNgType::Readable&           io,
               LanczosPlusPlus::LanczosOptions& lanczosOptions);

template <typename ModelType>
void mainLoop(InputNgType::Readable&           io,
              const ModelType&                 model,
              LanczosPlusPlus::LanczosOptions& lanczosOptions);

template <typename ComplexOrRealType>
void mainLoop0(InputNgType::Readable& io, LanczosPlusPlus::LanczosOptions& lanczosOptions);

#endif // LANCZOSDRIVER_H
