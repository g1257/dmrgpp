#ifndef MANYOMEGAS_H
#define MANYOMEGAS_H
#include "../../dmrgpp/src/Engine/InputCheck.h"
#include "../../dmrgpp/src/Engine/Qn.h"
#include "../../dmrgpp/src/Engine/ProgramGlobals.h"
#include "../../dmrgpp/src/Engine/SuperGeometry.h"
#include "../../dmrgpp/src/Engine/ParametersDmrgSolver.h"
#include "../../dmrgpp/src/Engine/ModelSelector.h"
#include "../../dmrgpp/src/Engine/ModelHelperLocal.h"
#include "../../dmrgpp/src/Engine/MatrixVectorKron/MatrixVectorKron.h"
#include "../../dmrgpp/src/Engine/MatrixVectorOnTheFly.h"
#include "../../dmrgpp/src/Engine/MatrixVectorStored.h"
#include "../../dmrgpp/src/Engine/LeftRightSuper.h"
#include "../../dmrgpp/src/Engine/BasisWithOperators.h"
#include "../../dmrgpp/src/Engine/DmrgSolver.h"
#include "../../dmrgpp/src/Engine/VectorWithOffset.h"

#include "PsiBase64.h"
#include "InputNg.h"
#include "LanczosSolver.h"
#include "Vector.h"
#include "ParamsDmftSolver.h"
namespace Dmft {

template<typename ComplexOrRealType>
class ManyOmegas {

public:

	typedef typename PsimagLite::Real<ComplexOrRealType>::Type RealType;
	typedef typename PsimagLite::Vector<RealType>::Type VectorRealType;
	typedef PsimagLite::InputNg<Dmrg::InputCheck> InputNgType;
	typedef Dmrg::ParametersDmrgSolver<RealType, InputNgType::Readable, Dmrg::Qn>
	ParametersDmrgSolverType;
	typedef Dmrg::SuperGeometry<ComplexOrRealType,
	        InputNgType::Readable,
	        Dmrg::ProgramGlobals> SuperGeometryType;
	typedef Dmrg::VectorWithOffset<ComplexOrRealType, Dmrg::Qn> VectorWithOffsetType;

	ManyOmegas manyOmegas(PsimagLite::String inputFile,
	                      RealType precision,
	                      bool echoInput)
	    : inputfile_(inputFile), precision_(precision), echoInput_(echoInput)
	{
		InputNgType::Writeable::readFile(data_, inputFile);
	}

	void run(bool dryRun)
	{
		for (SizeType i = offset; i < total; ++i) {
			RealType omega = omegas_[i];
			PsimagLite::String data2 = modifyOmega(data_, omega);
			PsimagLite::String sOptions = "";
			doOneRun(data2, sOptions);
		}
	}

	PsimagLite::String inputFile_;
	RealType precision_;
	bool echoInput_;
	PsimagLite::String data_;
};
}
#endif // MANYOMEGAS_H
