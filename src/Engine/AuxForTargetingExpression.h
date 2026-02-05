#ifndef AUXFORTARGETINGEXPRESSION_H
#define AUXFORTARGETINGEXPRESSION_H
#include "GetBraOrKet.h"
#include "GroupOfOneTimeEvolutions.h"
#include "InputCheck.h"
#include "InputNg.h"
#include "ProgramGlobals.h"
#include "Pvectors.h"
#include "Vector.h"

namespace Dmrg {

template <typename TargetingBaseType> class AuxForTargetingExpression {

public:

	using VectorWithOffsetType        = typename TargetingBaseType::VectorWithOffsetType;
	using ModelType                   = typename TargetingBaseType::ModelType;
	using ApplyOperatorExpressionType = typename TargetingBaseType::ApplyOperatorExpressionType;
	using VectorVectorWithOffsetType  = typename PsimagLite::Vector<VectorWithOffsetType>::Type;
	typedef
	    typename PsimagLite::Vector<typename PsimagLite::Vector<VectorWithOffsetType*>::Type>::
	        Type VectorVectorVectorWithOffsetType;
	using LeftRightSuperType           = typename ModelType::LeftRightSuperType;
	using VectorStringType             = PsimagLite::Vector<PsimagLite::String>::Type;
	using TargetParamsType             = typename TargetingBaseType::TargetParamsType;
	using InputValidatorType           = PsimagLite::InputNg<InputCheck>::Readable;
	using ComplexOrRealType            = typename VectorWithOffsetType::value_type;
	using RealType                     = typename PsimagLite::Real<ComplexOrRealType>::Type;
	using PvectorsType                 = Pvectors<TargetingBaseType>;
	using GroupOfOneTimeEvolutionsType = GroupOfOneTimeEvolutions<PvectorsType>;

	AuxForTargetingExpression(PvectorsType&                 pVectors,
	                          GroupOfOneTimeEvolutionsType& timeEvolve,
	                          ProgramGlobals::DirectionEnum dir,
	                          RealType                      Eg,
	                          SizeType                      currentCoo)
	    : pVectors_(pVectors)
	    , timeEvolve_(timeEvolve)
	    , direction_(dir)
	    , Eg_(Eg)
	    , currentCoo_(currentCoo)
	    , pIndexOutput_(0)
	{ }

	PvectorsType& pVectors() const { return pVectors_; }

	GroupOfOneTimeEvolutionsType& timeEvolve() const { return timeEvolve_; }

	void setPindexOutput(SizeType x) { pIndexOutput_ = x; }

	ProgramGlobals::DirectionEnum direction() const { return direction_; }

	const RealType& Eg() const { return Eg_; }

	SizeType currentCoO() const { return currentCoo_; }

	const SizeType pIndexOutput() const { return pIndexOutput_; }

	VectorWithOffsetType& getCurrentVectorNonConst(PsimagLite::String braOrKet) const
	{
		PsimagLite::GetBraOrKet getBraOrKet(braOrKet);
		if (getBraOrKet.isRvector()) {
			const SizeType pIndex = getBraOrKet.pIndex();
			if (pIndex >= tempVectors_.size())
				err("getCurrentVectorNonConst: out of range for " + braOrKet
				    + "\n");
			return tempVectors_[pIndex];
		}

		throw PsimagLite::RuntimeError(
		    "getCurrentVectorNonConst: psi or tvs cannot be modified\n");
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

	PvectorsType&                      pVectors_;
	GroupOfOneTimeEvolutionsType&      timeEvolve_;
	ProgramGlobals::DirectionEnum      direction_;
	RealType                           Eg_;
	SizeType                           currentCoo_;
	SizeType                           pIndexOutput_;
	mutable VectorVectorWithOffsetType tempVectors_;
	mutable VectorStringType           tempNames_;
};

}
#endif // AUXFORTARGETINGEXPRESSION_H
