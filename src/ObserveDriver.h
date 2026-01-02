#ifndef OBSERVEDRIVER_H
#define OBSERVEDRIVER_H

#include <unistd.h>
#define USE_PTHREADS_OR_NOT_NG
#include "BasisWithOperators.h"
#include "CrsMatrix.h"
#include "DmrgSolver.h" // only used for types
#include "Geometry/Geometry.h"
#include "InputCheck.h"
#include "InputFromDataOrNot.h"
#include "InputNg.h"
#include "Io/IoSelector.h"
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
#include "VectorWithOffset.h"
#include "VectorWithOffsets.h"

namespace Dmrg {

typedef PsimagLite::IoSelector::In IoInputType;

#ifndef USE_FLOAT
typedef double RealType;
#else
typedef float RealType;
#endif

typedef std::complex<RealType> ComplexType;

typedef PsimagLite::CrsMatrix<ComplexType> MySparseMatrixComplex;
typedef PsimagLite::CrsMatrix<RealType>    MySparseMatrixReal;

typedef PsimagLite::InputNg<InputCheck>                                 InputNgType;
typedef ParametersDmrgSolver<RealType, InputNgType::Readable, Dmrg::Qn> ParametersDmrgSolverType;

template <typename VectorWithOffsetType, typename ModelType>
bool observeOneFullSweep(IoInputType&              io,
                         const ModelType&          model,
                         const PsimagLite::String& list,
                         SizeType                  orbitals);
}

#endif // OBSERVEDRIVER_H
