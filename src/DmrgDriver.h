#ifndef DMRGDRIVER_H
#define DMRGDRIVER_H
#define USE_PTHREADS_OR_NOT_NG
#include "DmrgSolver.h"
#include "InputNg.h"
#include "InputCheck.h"
#include "ParametersDmrgSolver.h"
#include "ModelSelector.h"
#include "Geometry/Geometry.h"
#include "ModelHelperLocal.h"
#include "ModelHelperSu2.h"
#include "MatrixVectorOnTheFly.h"
#include "MatrixVectorStored.h"
#include "MatrixVectorKron/MatrixVectorKron.h"
#include "TargetingBase.h"
#include "VectorWithOffset.h"
#include "VectorWithOffsets.h"
#include "BasisWithOperators.h"
#include "LeftRightSuper.h"
#include "LanczosSolver.h"
#include "ChebyshevSolver.h"
#include "Operators.h"
#include "CrsMatrix.h"
#include "OperatorSpec.h"
#include "CanonicalExpression.h"

#ifndef USE_FLOAT
typedef double RealType;
#else
typedef float RealType;
#endif

struct OperatorOptions {

	OperatorOptions()
	    : site(0),
	      dof(0),
	      label(""),
	      opexpr(""),
	      hasOperatorExpression(false),
	      transpose(false),
	      enabled(false)
	{}

	SizeType site;
	SizeType dof;
	PsimagLite::String label;
	PsimagLite::String opexpr;
	bool hasOperatorExpression;
	bool transpose;
	bool enabled;
};

typedef PsimagLite::InputNg<Dmrg::InputCheck> InputNgType;
typedef Dmrg::ParametersDmrgSolver<RealType,InputNgType::Readable> ParametersDmrgSolverType;

void usageOperator();

template<typename ModelBaseType>
void operatorDriver(const ModelBaseType& model, const OperatorOptions& obsOptions)
{
	typedef typename ModelBaseType::ModelHelperType ModelHelperType;
	typedef typename ModelHelperType::OperatorsType OperatorsType;
	typedef typename OperatorsType::OperatorType OperatorType;
	typedef Dmrg::OperatorSpec<ModelBaseType> OperatorSpecType;

	if (obsOptions.hasOperatorExpression && obsOptions.label != "") {
		std::cerr<<"You must provide exactly one option: -l or -e;";
		std::cerr<<" both were given\n";
		usageOperator();
		return;
	}

	if (!obsOptions.hasOperatorExpression && obsOptions.label == "") {
		std::cerr<<"You must provide exactly one option: -l or -e;";
		std::cerr<<" none were given\n";
		usageOperator();
		return;
	}

	OperatorType opC;

	if (obsOptions.hasOperatorExpression) {
		OperatorSpecType opSpec(model);
		int site = -1;
		PsimagLite::CanonicalExpression<OperatorSpecType> canonicalExpression(opSpec);
		opC = canonicalExpression(obsOptions.opexpr, site);
	} else {
		if (obsOptions.label == "B") {
			model.printBasis(obsOptions.site);
			return;
		}

		opC = model.naturalOperator(obsOptions.label,
		                            obsOptions.site,
		                            obsOptions.dof);
		std::cerr<<"label="<<obsOptions.label<<" site="<<obsOptions.site;
		std::cerr<<" dof="<<obsOptions.dof<<"\n";
	}

	if (obsOptions.transpose) opC.dagger();

	opC.write(std::cout);
}

template<typename SolverType, typename VectorWithOffsetType>
void mainLoop4(typename SolverType::LanczosMatrixType::ModelType::GeometryType&,
               const ParametersDmrgSolverType&,
               InputNgType::Readable&,
               const OperatorOptions&,
               PsimagLite::String);

#endif // DMRGDRIVER_H

