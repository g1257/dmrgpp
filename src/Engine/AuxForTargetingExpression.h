#ifndef AUXFORTARGETINGEXPRESSION_H
#define AUXFORTARGETINGEXPRESSION_H
#include "Vector.h"
#include "GetBraOrKet.h"
#include "ProgramGlobals.h"

namespace Dmrg {

template<typename TargetingBaseType>
class AuxForTargetingExpression {

public:

	typedef typename TargetingBaseType::VectorWithOffsetType VectorWithOffsetType;
	typedef typename TargetingBaseType::ModelType ModelType;
	typedef typename TargetingBaseType::ApplyOperatorExpressionType ApplyOperatorExpressionType;
	typedef typename PsimagLite::Vector<VectorWithOffsetType>::Type VectorVectorWithOffsetType;
	typedef typename PsimagLite::Vector<typename
	PsimagLite::Vector<VectorWithOffsetType*>::Type>::Type VectorVectorVectorWithOffsetType;
	typedef typename ModelType::LeftRightSuperType LeftRightSuperType;

	AuxForTargetingExpression(const ApplyOperatorExpressionType& aoe,
	                          const ModelType& model,
	                          const LeftRightSuperType& lrs,
	                          ProgramGlobals::DirectionEnum dir)
	    : aoe_(aoe), model_(model), lrs_(lrs), direction_(dir)
	{}

	const ApplyOperatorExpressionType& aoe() const { return aoe_; }

	const ModelType& model() const { return model_; }

	const LeftRightSuperType& lrs() const { return lrs_; }

	const ProgramGlobals::DirectionEnum direction() const { return direction_; }

	VectorWithOffsetType& getCurrentVector(PsimagLite::String braOrKet) const
	{
		throw PsimagLite::RuntimeError("getCurrentVector\n");
//		PsimagLite::GetBraOrKet getBraOrKet(braOrKet);

//		if (getBraOrKet.isPvector()) {
//			const SizeType pIndex = getBraOrKet.pIndex();
//			if (pIndex >= pvectors.size())
//				err("getVector: out of range for " + braOrKet + "\n");
//			return pvectors[pIndex];
//		}

//		const SizeType sectorIndex = getBraOrKet.sectorIndex();
//		return *(psi[sectorIndex][getBraOrKet.levelIndex()]);
	}

	void createTemporaryVector(PsimagLite::String str) const
	{
		err("createTemporaryVector\n");
	}

private:

	const ApplyOperatorExpressionType& aoe_;
	const ModelType& model_;
	const LeftRightSuperType lrs_;
	ProgramGlobals::DirectionEnum direction_;
};

}
#endif // AUXFORTARGETINGEXPRESSION_H
