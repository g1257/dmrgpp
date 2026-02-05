#ifndef TARGETSELECTOR_H
#define TARGETSELECTOR_H
#include "TargetingBase.h"
#include "TargetingGroundState.h"
// start headers: // DO NOT REMOVE MARK
#include "TargetingCVEvolution.h"
#include "TargetingChebyshev.h"
#include "TargetingCorrection.h"
#include "TargetingCorrectionVector.h"
#include "TargetingDynamic.h"
#include "TargetingExpression.h"
#include "TargetingMetts.h"
#include "TargetingRixsDynamic.h"
#include "TargetingRixsStatic.h"
#include "TargetingTimeStep.h"
// end targets DO NOT REMOVE MARK

namespace Dmrg {

template <typename TargetingBaseType> class TargetSelector {

	using MatrixVectorType       = typename TargetingBaseType::MatrixVectorType;
	using CheckpointType         = typename TargetingBaseType::CheckpointType;
	using ModelType              = typename MatrixVectorType::ModelType;
	using ParametersType         = typename ModelType::ParametersType;
	using RealType               = typename ModelType::RealType;
	using ModelHelperType        = typename ModelType::ModelHelperType;
	using LeftRightSuperType     = typename ModelHelperType::LeftRightSuperType;
	using WaveFunctionTransfType = typename TargetingBaseType::WaveFunctionTransfType;
	using BasisWithOperatorsType = typename LeftRightSuperType::BasisWithOperatorsType;
	using BasisType              = typename BasisWithOperatorsType::BasisType;
	using QnType                 = typename BasisType::QnType;
	using OptionsType            = typename ParametersType::OptionsType;
	using InputValidatorType     = typename ModelType::InputValidatorType;
	using LanczosSolverType      = typename TargetingBaseType::LanczosSolverType;
	using VectorWithOffsetType   = typename TargetingBaseType::VectorWithOffsetType;
	using TargetingGroundStateType
	    = TargetingGroundState<LanczosSolverType, VectorWithOffsetType>;
	using VectorStringType = PsimagLite::Vector<PsimagLite::String>::Type;
	// start targets here:  DO NOT REMOVE MARK
	using TargetingTimeStepType  = TargetingTimeStep<LanczosSolverType, VectorWithOffsetType>;
	using TargetingChebyshevType = TargetingChebyshev<LanczosSolverType, VectorWithOffsetType>;
	using TargetingDynamicType   = TargetingDynamic<LanczosSolverType, VectorWithOffsetType>;
	using TargetingCorrectionVectorType
	    = TargetingCorrectionVector<LanczosSolverType, VectorWithOffsetType>;
	using TargetingCorrectionType
	    = TargetingCorrection<LanczosSolverType, VectorWithOffsetType>;
	using TargetingMettsType = TargetingMetts<LanczosSolverType, VectorWithOffsetType>;
	using TargetingRixsStaticType
	    = TargetingRixsStatic<LanczosSolverType, VectorWithOffsetType>;
	using TargetingRixsDynamicType
	    = TargetingRixsDynamic<LanczosSolverType, VectorWithOffsetType>;
	using TargetingExpressionType
	    = TargetingExpression<LanczosSolverType, VectorWithOffsetType>;
	using TargetingCVEvolutionType
	    = TargetingCVEvolution<LanczosSolverType, VectorWithOffsetType>;
	// end targets DO NOT REMOVE MARK

public:

	TargetSelector(const LeftRightSuperType&            lrs,
	               const CheckpointType&                checkPoint,
	               const WaveFunctionTransfType&        wft,
	               const typename QnType::VectorQnType& quantumSector,
	               InputValidatorType&                  ioIn)
	    : psi_(nullptr)
	    , lrs_(lrs)
	    , checkPoint_(checkPoint)
	    , wft_(wft)
	    , quantumSector_(quantumSector)
	    , ioIn_(ioIn)
	{ }

	~TargetSelector()
	{
		delete psi_;
		psi_ = nullptr;
	}

	TargetSelector(const TargetSelector&) = delete;

	TargetSelector& operator=(const TargetSelector&) = delete;

	TargetingBaseType& operator()()
	{
		if (psi_)
			return *psi_;

		PsimagLite::String targeting = getTargeting(checkPoint_.model().params().options);
		return operator()(targeting);
	}

	TargetingBaseType& operator()(PsimagLite::String targeting)
	{
		if (psi_)
			err("TargetingBaseType::operator(): can only be called multiple times"
			    + PsimagLite::String(" without argument\n"));

		check1(targeting);

		assert(0 < quantumSector_.size());
		const QnType& qn = quantumSector_[0];
		if (targeting == "GroundStateTargeting") {
			psi_ = new TargetingGroundStateType(lrs_, checkPoint_, wft_, qn, ioIn_);
			// named targets start  DO NOT REMOVE MARK
		} else if (targeting == "TimeStepTargeting" || targeting == "TargetingAncilla") {
			psi_ = new TargetingTimeStepType(
			    lrs_, checkPoint_, wft_, qn, ioIn_, targeting);
		} else if (targeting == "DynamicTargeting") {
			psi_ = new TargetingDynamicType(lrs_, checkPoint_, wft_, qn, ioIn_);
		} else if (targeting == "CorrectionVectorTargeting") {
			psi_
			    = new TargetingCorrectionVectorType(lrs_, checkPoint_, wft_, qn, ioIn_);
		} else if (targeting == "TargetingChebyshev") {
			psi_ = new TargetingChebyshevType(lrs_, checkPoint_, wft_, qn, ioIn_);
		} else if (targeting == "CorrectionTargeting") {
			psi_ = new TargetingCorrectionType(lrs_, checkPoint_, wft_, qn, ioIn_);
		} else if (targeting == "MettsTargeting") {
			psi_ = new TargetingMettsType(lrs_, checkPoint_, wft_, qn, ioIn_);
		} else if (targeting == "TargetingRixsStatic") {
			psi_ = new TargetingRixsStaticType(lrs_, checkPoint_, wft_, qn, ioIn_);
		} else if (targeting == "TargetingRixsDynamic") {
			psi_ = new TargetingRixsDynamicType(lrs_, checkPoint_, wft_, qn, ioIn_);
		} else if (targeting == "TargetingExpression") {
			psi_ = new TargetingExpressionType(lrs_, checkPoint_, wft_, qn, ioIn_);
		} else if (targeting == "TargetingCVEvolution") {
			psi_ = new TargetingCVEvolutionType(lrs_, checkPoint_, wft_, qn, ioIn_);
			// end targets DO NOT REMOVE MARK
		} else {
			err("Unknown targeting " + targeting + "\n");
		}

		psi_->postCtor();

		return *psi_;
	}

private:

	PsimagLite::String getTargeting(const OptionsType& options) const
	{
		PsimagLite::String targeting("GroundStateTargeting");

		VectorStringType targets = { "GroundStateTargeting",
			                     "TimeStepTargeting",
			                     "AdaptiveDynamicTargeting",
			                     "DynamicTargeting",
			                     "CorrectionVectorTargeting",
			                     "CorrectionTargeting",
			                     "MettsTargeting",
			                     "TargetingAncilla",
			                     "TargetingCorrelations",
			                     "TargetingInSitu",
			                     "TargetingRixsStatic",
			                     "TargetingRixsDynamic",
			                     "TargetingChebyshev",
			                     "TargetingExpression",
			                     "TargetingMultiQ",
			                     "TargetingCVEvolution" };

		const SizeType totalTargets = targets.size();
		SizeType       count        = 0;
		for (SizeType i = 0; i < totalTargets; ++i) {
			if (options.isSet(targets[i])) {
				if (targeting == "AdaptiveDynamicTargeting"
				    && targets[i] == "DynamicTargeting")
					continue;
				targeting = targets[i];
				count++;
			}
		}

		if (count == 1)
			return targeting;

		if (count > 1)
			err("Only one targeting at a time supported\n");

		std::cerr << " No explicit targeting found, asumming " << targeting << "\n";
		std::cout << " No explicit targeting found, asumming " << targeting << "\n";

		return targeting;
	}

	void check1(PsimagLite::String) const { }

	TargetingBaseType*                   psi_;
	const LeftRightSuperType&            lrs_;
	const CheckpointType&                checkPoint_;
	const WaveFunctionTransfType&        wft_;
	const typename QnType::VectorQnType& quantumSector_;
	InputValidatorType&                  ioIn_;
};
}
#endif // TARGETSELECTOR_H
