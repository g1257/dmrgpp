#ifndef AUXFORTARGETINGEXPRESSION_H
#define AUXFORTARGETINGEXPRESSION_H
#include "Vector.h"
#include "GetBraOrKet.h"
#include "ProgramGlobals.h"
#include "InputNg.h"
#include "InputCheck.h"
#include "Pvectors.h"
#include "TimeEvolveForTargetingExpression.h"

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
	typedef Pvectors<TargetingBaseType> PvectorsType;
	typedef TimeEvolveForTargetingExpression<PvectorsType> TimeEvolveForTargetingExpressionType;

	AuxForTargetingExpression(PvectorsType& pVectors,
	                          TimeEvolveForTargetingExpressionType& timeEvolve,
	                          ProgramGlobals::DirectionEnum dir,
	                          RealType Eg)
	    : pVectors_(pVectors),
	      timeEvolve_(timeEvolve),
	      direction_(dir),
	      Eg_(Eg),
	      pIndexOutput_(0)
	{}

	PvectorsType& pVectors() const { return pVectors_; }

	TimeEvolveForTargetingExpressionType& timeEvolve() const { return timeEvolve_; }

	void setPindexOutput(SizeType x) { pIndexOutput_ = x; }

	ProgramGlobals::DirectionEnum direction() const { return direction_; }

	const RealType& Eg() const { return Eg_; }

	const SizeType pIndexOutput() const { return pIndexOutput_; }

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

	const VectorVectorWithOffsetType& tempVectors() const { return tempVectors_; }

	const VectorStringType& tempNames() const { return tempNames_; }

private:

	PvectorsType& pVectors_;
	TimeEvolveForTargetingExpressionType& timeEvolve_;
	ProgramGlobals::DirectionEnum direction_;
	RealType Eg_;
	SizeType pIndexOutput_;
	mutable VectorVectorWithOffsetType tempVectors_;
	mutable VectorStringType tempNames_;
};

}
#endif // AUXFORTARGETINGEXPRESSION_H
