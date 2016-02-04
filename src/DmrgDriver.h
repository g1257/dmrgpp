#ifndef DMRGDRIVER_H
#define DMRGDRIVER_H
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
#include "MatrixVectorKron.h"
#include "TargetingBase.h"
#include "MettsTargetting.h"
#include "VectorWithOffset.h"
#include "VectorWithOffsets.h"
#include "BasisWithOperators.h"
#include "LeftRightSuper.h"
#include "LanczosSolver.h"
#include "ChebyshevSolver.h"
#include "Operators.h"
#include "CrsMatrix.h"

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
	      fermionicSign(0),
	      transpose(false),
	      enabled(false)
	{}

	SizeType site;
	SizeType dof;
	PsimagLite::String label;
	int fermionicSign;
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
	typedef typename ModelHelperType::SparseElementType SparseElementType;
	typedef PsimagLite::Matrix<SparseElementType> MatrixType;
	typedef std::pair<SizeType,SizeType> PairType;

	if (obsOptions.label=="" || obsOptions.fermionicSign == 0) {
		usageOperator();
		return;
	}

	MatrixType opC = model.naturalOperator(obsOptions.label,
	                                       obsOptions.site,
	                                       obsOptions.dof);
	std::cerr<<"#label="<<obsOptions.label<<" site="<<obsOptions.site;
	std::cerr<<" dof="<<obsOptions.dof<<"\n";

	Dmrg::Su2Related su2Related;

	MatrixType opC2;
	if (obsOptions.transpose)
		transposeConjugate(opC2,opC);

	OperatorType opC3((obsOptions.transpose) ? opC2 : opC,
	                  obsOptions.fermionicSign,
	                  PairType(0,0),
	                  1,
	                  su2Related);

	opC3.save(std::cout);
}

template<typename SolverType, typename VectorWithOffsetType>
void mainLoop4(typename SolverType::LanczosMatrixType::ModelType::GeometryType&,
               const ParametersDmrgSolverType&,
               InputNgType::Readable&,
               const OperatorOptions&,
               PsimagLite::String);

#endif // DMRGDRIVER_H

