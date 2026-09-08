#ifndef OBSERVEDRIVER_H
#define OBSERVEDRIVER_H

#include "BasisWithOperators.h"
#include "DmrgSolver.h" // only used for types
#include "InputCheck.h"
#include "InputFromDataOrNot.h"
#include "LeftRightSuper.h"
#include "ModelBase.h"
#include "ModelHelperLocal.h"
#include "ModelSelector.h"
#include "ObservableLibrary.h"
#include "Observer.h"
#include "Operators.h"
#include "Provenance.h"
#include "SuperGeometry.h"
#include "TargetingCorrection.h"
#include "TargetingCorrectionVector.h"
#include "TargetingDynamic.h"
#include "TargetingGroundState.h"
#include "TargetingMetts.h"
#include "TargetingTimeStep.h"
#include "VectorWithOffsets.h"
#include <PsimagLite/CrsMatrix.h>
#include <PsimagLite/Geometry/Geometry.h>
#include <PsimagLite/InputNg.h>
#include <PsimagLite/Io/IoSelector.h>
#include <unistd.h>

namespace Dmrg {

using IoInputType = PsimagLite::IoSelector::In;

#ifndef USE_FLOAT
using RealType = double;
#else
using RealType = float;
#endif

using ComplexType = std::complex<RealType>;

using MySparseMatrixComplex = PsimagLite::CrsMatrix<ComplexType>;
using MySparseMatrixReal    = PsimagLite::CrsMatrix<RealType>;

using InputNgType              = PsimagLite::InputNg<InputCheck>;
using ParametersDmrgSolverType = ParametersDmrgSolver<RealType, InputNgType::Readable, Dmrg::Qn>;

template <typename VectorWithOffsetType, typename ModelType>
bool observeOneFullSweep(IoInputType&              io,
                         const ModelType&          model,
                         const PsimagLite::String& list,
                         SizeType                  orbitals);
}

#endif // OBSERVEDRIVER_H
