#ifndef AUXFORTARGETINGEXPRESSION_H
#define AUXFORTARGETINGEXPRESSION_H

namespace Dmrg {

template<typename TargetingBaseType>
struct AuxForTargetingExpression {

	typedef typename TargetingBaseType::VectorWithOffsetType VectorWithOffsetType;
	typedef typename TargetingBaseType::ModelType ModelType;
	typedef typename TargetingBaseType::ApplyOperatorExpressionType ApplyOperatorExpressionType;
	typedef typename PsimagLite::Vector<VectorWithOffsetType>::Type VectorVectorWithOffsetType;
	typedef typename PsimagLite::Vector<typename
	PsimagLite::Vector<VectorWithOffsetType*>::Type>::Type VectorVectorVectorWithOffsetType;
	typedef typename ModelType::LeftRightSuperType LeftRightSuperType;

	AuxForTargetingExpression(const ApplyOperatorExpressionType& aoe_,
	                          const ModelType& model_,
	                          const LeftRightSuperType& lrs_,
	                          const VectorVectorVectorWithOffsetType& psi_,
	                          const VectorVectorWithOffsetType& pvectors_,
	                          ProgramGlobals::DirectionEnum dir)
	    : aoe(aoe_), model(model_), lrs(lrs_), psi(psi_), pvectors(pvectors_), direction(dir)
	{}

	const ApplyOperatorExpressionType& aoe;
	const ModelType& model;
	const LeftRightSuperType lrs;
	const VectorVectorVectorWithOffsetType& psi;
	const VectorVectorWithOffsetType& pvectors;
	ProgramGlobals::DirectionEnum direction;
};

}
#endif // AUXFORTARGETINGEXPRESSION_H
