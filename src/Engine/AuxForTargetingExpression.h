#ifndef AUXFORTARGETINGEXPRESSION_H
#define AUXFORTARGETINGEXPRESSION_H
#include "Vector.h"
#include "GetBraOrKet.h"
#include "ProgramGlobals.h"
#include "InputNg.h"
#include "InputCheck.h"

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
	typedef PsimagLite::Vector<PsimagLite::String>::Type VectorStringType;
	typedef typename TargetingBaseType::TargetParamsType TargetParamsType;
	typedef PsimagLite::InputNg<InputCheck>::Readable InputValidatorType;
	typedef typename VectorWithOffsetType::value_type ComplexOrRealType;
	typedef typename PsimagLite::Real<ComplexOrRealType>::Type RealType;

	AuxForTargetingExpression(const ApplyOperatorExpressionType& aoe,
	                          const LeftRightSuperType& lrs,
	                          ProgramGlobals::DirectionEnum dir,
	                          const TargetParamsType& tstStruct,
	                          RealType Eg)
	    : aoe_(aoe), lrs_(lrs), direction_(dir), tstStruct_(tstStruct), Eg_(Eg)
	{}

	const ApplyOperatorExpressionType& aoe() const { return aoe_; }

	ApplyOperatorExpressionType& aoeNonConst()
	{
		ApplyOperatorExpressionType* aoePtr = const_cast<ApplyOperatorExpressionType*>(&aoe_);
		return *aoePtr;
	}

	const LeftRightSuperType& lrs() const { return lrs_; }

	ProgramGlobals::DirectionEnum direction() const { return direction_; }

	const TargetParamsType& tstStruct() const { return tstStruct_; }

	const RealType& Eg() const { return Eg_; }

	const VectorWithOffsetType& getCurrentVectorConst(PsimagLite::String braOrKet) const
	{
		PsimagLite::GetBraOrKet getBraOrKet(braOrKet);
		if (getBraOrKet.isPvector()) {
			const SizeType pIndex = getBraOrKet.pIndex();
			if (pIndex >= aoe_.targetVectors().size())
				err("getVector: out of range for " + braOrKet + "\n");
			return aoe_.targetVectors()[pIndex];
		} else if (getBraOrKet.isRvector()) {
			throw PsimagLite::RuntimeError("reserved vector\n");
		}

		const SizeType sectorIndex = getBraOrKet.sectorIndex();
		return *(aoe_.psiConst()[sectorIndex][getBraOrKet.levelIndex()]);
	}

	VectorWithOffsetType& getCurrentVectorNonConst(PsimagLite::String braOrKet) const
	{
		PsimagLite::GetBraOrKet getBraOrKet(braOrKet);
		if (getBraOrKet.isRvector()) {
			const SizeType pIndex = getBraOrKet.pIndex();
			if (pIndex >= tempVectors_.size())
				err("getCurrentVectorNonConst: out of range for " + braOrKet + "\n");
			return tempVectors_[pIndex];
		}

		throw PsimagLite::RuntimeError("getCurrentVectorNonConst: psi or tvs cannot be modified\n");
	}

	PsimagLite::String createTemporaryVector(PsimagLite::String str) const
	{
		const SizeType n = tempVectors_.size();
		tempVectors_.push_back(VectorWithOffsetType());
		tempNames_.push_back(str);
		return "R" + ttos(n);
	}

	const VectorVectorWithOffsetType& tempVectors() const
	{
		return tempVectors_;
	}

	const VectorStringType& tempNames() const
	{
		return tempNames_;
	}

private:

	const ApplyOperatorExpressionType& aoe_;
	const LeftRightSuperType lrs_;
	ProgramGlobals::DirectionEnum direction_;
	const TargetParamsType& tstStruct_;
	RealType Eg_;
	mutable VectorVectorWithOffsetType tempVectors_;
	mutable VectorStringType tempNames_;
};

}
#endif // AUXFORTARGETINGEXPRESSION_H
