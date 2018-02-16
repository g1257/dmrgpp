#ifndef OBSERVEDRIVER_H
#define OBSERVEDRIVER_H

#include <unistd.h>
#define USE_PTHREADS_OR_NOT_NG
#include "Observer.h"
#include "ObservableLibrary.h"
#include "IoSelector.h"
#include "Operators.h"
#include "Geometry/Geometry.h"
#include "CrsMatrix.h"
#include "ModelHelperLocal.h"
#include "ModelHelperSu2.h"
#include "VectorWithOffset.h"
#include "VectorWithOffsets.h"
#include "DmrgSolver.h" // only used for types
#include "TargetingGroundState.h"
#include "TargetingTimeStep.h"
#include "TargetingDynamic.h"
#include "TargetingAdaptiveDynamic.h"
#include "TargetingCorrection.h"
#include "TargetingCorrectionVector.h"
#include "TargetingMetts.h"
#include "BasisWithOperators.h"
#include "LeftRightSuper.h"
#include "InputNg.h"
#include "Provenance.h"
#include "InputCheck.h"
#include "ModelSelector.h"
#include "ArchiveFiles.h"
#include "CvectorSize.h"

namespace Dmrg {

typedef PsimagLite::IoSelector::In IoInputType;

#ifndef USE_FLOAT
typedef double RealType;
#else
typedef float RealType;
#endif

typedef std::complex<RealType> ComplexType;

typedef  PsimagLite::CrsMatrix<ComplexType> MySparseMatrixComplex;
typedef  PsimagLite::CrsMatrix<RealType> MySparseMatrixReal;

typedef PsimagLite::InputNg<InputCheck> InputNgType;
typedef ParametersDmrgSolver<RealType,InputNgType::Readable> ParametersDmrgSolverType;

template<typename VectorWithOffsetType,
         typename ModelType>
bool observeOneFullSweep(IoInputType& io,
                         const ModelType& model,
                         const PsimagLite::String& list,
                         bool hasTimeEvolution,
                         SizeType orbitals);
}

#endif // OBSERVEDRIVER_H
