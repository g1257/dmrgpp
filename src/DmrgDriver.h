#ifndef DMRGDRIVER_H
#define DMRGDRIVER_H
#define USE_PTHREADS_OR_NOT_NG
#include "BasisWithOperators.h"
#include "CanonicalExpression.h"
#include "ChebyshevSolver.h"
#include "CrsMatrix.h"
#include "DmrgSolver.h"
#include "InputCheck.h"
#include "InputNg.h"
#include "LanczosSolver.h"
#include "LeftRightSuper.h"
#include "MatrixVectorKron/MatrixVectorKron.h"
#include "MatrixVectorOnTheFly.h"
#include "MatrixVectorStored.h"
#include "ModelBase.h"
#include "ModelHelperLocal.h"
#include "ModelSelector.h"
#include "OperatorSpec.h"
#include "Operators.h"
#include "ParametersDmrgSolver.h"
#include "ProgramGlobals.h"
#include "SuperGeometry.h"
#include "TargetingBase.h"
#include "VectorWithOffset.h"
#include "VectorWithOffsets.h"

#ifndef USE_FLOAT
typedef double RealType;
#else
typedef float RealType;
#endif

struct OperatorOptions {

	OperatorOptions()
	    : site(0)
	    , dof(0)
	    , label("")
	    , opexpr("")
	    , hasOperatorExpression(false)
	    , transpose(false)
	    , enabled(false)
	{
	}

	SizeType site;
	SizeType dof;
	PsimagLite::String label;
	PsimagLite::String opexpr;
	bool hasOperatorExpression;
	bool transpose;
	bool enabled;
};

typedef PsimagLite::InputNg<Dmrg::InputCheck> InputNgType;
typedef Dmrg::ParametersDmrgSolver<RealType,
    InputNgType::Readable,
    Dmrg::Qn>
    ParametersDmrgSolverType;

void usageOperator();

template <typename ModelBaseType>
void operatorDriver(const ModelBaseType& model, const OperatorOptions& obsOptions)
{
	typedef typename ModelBaseType::ModelHelperType ModelHelperType;
	typedef typename ModelHelperType::OperatorsType OperatorsType;
	typedef typename OperatorsType::OperatorType OperatorType;
	typedef Dmrg::OperatorSpec<ModelBaseType, OperatorType> OperatorSpecType;

	if (obsOptions.hasOperatorExpression && obsOptions.label != "") {
		std::cerr << "You must provide exactly one option: -l or -e;";
		std::cerr << " both were given\n";
		usageOperator();
		return;
	}

	if (!obsOptions.hasOperatorExpression && obsOptions.label == "") {
		if (model.introspect())
			return;

		std::cerr << "You must provide exactly one option: -l or -e;";
		std::cerr << " none were given\n";
		usageOperator();
		return;
	}

	OperatorType opC;
	const OperatorType opEmpty;

	if (obsOptions.hasOperatorExpression) {
		OperatorSpecType opSpec(model);
		int site = -1;
		PsimagLite::CanonicalExpression<OperatorSpecType> canonicalExpression(opSpec);

		canonicalExpression(opC, obsOptions.opexpr, opEmpty, site);
	} else {
		if (obsOptions.label == "B") {
			model.printBasis(obsOptions.site);
			return;
		}

		if (obsOptions.label == "H") {
			model.printTerms();
			return;
		}

		opC = model.naturalOperator(obsOptions.label, obsOptions.site, obsOptions.dof);
		std::cerr << "label=" << obsOptions.label << " site=" << obsOptions.site;
		std::cerr << " dof=" << obsOptions.dof << "\n";
	}

	if (obsOptions.transpose)
		opC.dagger();

	opC.write(std::cout);
}

template <typename SolverType, typename VectorWithOffsetType>
void mainLoop4(typename SolverType::MatrixType::ModelType::SuperGeometryType&,
    const ParametersDmrgSolverType&,
    InputNgType::Readable&,
    const OperatorOptions&);

#endif // DMRGDRIVER_H
