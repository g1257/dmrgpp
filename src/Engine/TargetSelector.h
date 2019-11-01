#ifndef TARGETSELECTOR_H
#define TARGETSELECTOR_H
#include "TargetingBase.h"
#include "TargetingGroundState.h"
#include "TargetingTimeStep.h"
#include "TargetingDynamic.h"
#include "TargetingCorrection.h"
#include "TargetingCorrectionVector.h"
#include "TargetingChebyshev.h"
#include "TargetingMetts.h"
#include "TargetingRixsStatic.h"
#include "TargetingRixsDynamic.h"
#include "TargetingExpression.h"

namespace Dmrg {

template<typename TargetingBaseType>
class TargetSelector {

	typedef typename TargetingBaseType::MatrixVectorType MatrixVectorType;
	typedef typename MatrixVectorType::ModelType ModelType;
	typedef typename ModelType::ParametersType ParametersType;
	typedef typename ModelType::RealType RealType;
	typedef typename ModelType::ModelHelperType ModelHelperType;
	typedef typename ModelHelperType::LeftRightSuperType LeftRightSuperType;
	typedef typename TargetingBaseType::WaveFunctionTransfType WaveFunctionTransfType;
	typedef typename LeftRightSuperType::BasisWithOperatorsType BasisWithOperatorsType;
	typedef typename BasisWithOperatorsType::BasisType BasisType;
	typedef typename BasisType::QnType QnType;
	typedef typename ParametersType::OptionsType OptionsType;
	typedef typename ModelType::InputValidatorType InputValidatorType;
	typedef typename TargetingBaseType::LanczosSolverType LanczosSolverType;
	typedef typename TargetingBaseType::VectorWithOffsetType VectorWithOffsetType;
	typedef TargetingGroundState<LanczosSolverType,VectorWithOffsetType> TargetingGroundStateType;
	typedef TargetingTimeStep<LanczosSolverType,VectorWithOffsetType> TargetingTimeStepType;
	typedef TargetingChebyshev<LanczosSolverType,VectorWithOffsetType> TargetingChebyshevType;
	typedef TargetingDynamic<LanczosSolverType,VectorWithOffsetType> TargetingDynamicType;
	typedef TargetingCorrectionVector<LanczosSolverType,VectorWithOffsetType>
	TargetingCorrectionVectorType;
	typedef TargetingCorrection<LanczosSolverType,VectorWithOffsetType> TargetingCorrectionType;
	typedef TargetingMetts<LanczosSolverType,VectorWithOffsetType> TargetingMettsType;
	typedef TargetingRixsStatic<LanczosSolverType,VectorWithOffsetType> TargetingRixsStaticType;
	typedef TargetingRixsDynamic<LanczosSolverType,VectorWithOffsetType> TargetingRixsDynamicType;
	typedef TargetingExpression<LanczosSolverType,VectorWithOffsetType> TargetingExpressionType;
	typedef PsimagLite::Vector<PsimagLite::String>::Type VectorStringType;

public:

	TargetSelector(const LeftRightSuperType& lrs,
	               const ModelType& model,
	               const WaveFunctionTransfType& wft,
	               const typename QnType::VectorQnType& quantumSector,
	               InputValidatorType& ioIn)
	    : psi_(nullptr),
	      lrs_(lrs),
	      model_(model),
	      wft_(wft),
	      quantumSector_(quantumSector),
	      ioIn_(ioIn)
	{}

	~TargetSelector()
	{
		delete psi_;
		psi_ = nullptr;
	}

	TargetSelector(const TargetSelector&) = delete;

	TargetSelector& operator=(const TargetSelector&) = delete;

	TargetingBaseType& operator()()
	{
		if (psi_) return *psi_;

		PsimagLite::String targeting = getTargeting(model_.params().options);
		return operator()(targeting);
	}

	TargetingBaseType& operator()(PsimagLite::String targeting)
	{
		if (psi_)
			err("TargetingBaseType::operator(): can only be called multiple times" +
			    PsimagLite::String(" without argument\n"));

		check1(targeting);

		assert(0 < quantumSector_.size());
		const QnType& qn = quantumSector_[0];
		if (targeting=="TimeStepTargeting" || targeting == "TargetingAncilla") {
			psi_ = new TargetingTimeStepType(lrs_,model_,wft_,qn,ioIn_, targeting);
		} else if (targeting=="DynamicTargeting") {
			psi_ = new TargetingDynamicType(lrs_,model_,wft_,qn,ioIn_);
		} else if (targeting=="CorrectionVectorTargeting") {
			psi_ = new TargetingCorrectionVectorType(lrs_,model_,wft_,qn,ioIn_);
		} else if (targeting=="TargetingChebyshev") {
			psi_ = new TargetingChebyshevType(lrs_,model_,wft_,qn,ioIn_);
		} else if (targeting=="CorrectionTargeting") {
			psi_ = new TargetingCorrectionType(lrs_,model_,wft_,qn,ioIn_);
		} else if (targeting == "GroundStateTargeting") {
			psi_ = new TargetingGroundStateType(lrs_,model_,wft_,qn,ioIn_);
		} else if (targeting == "MettsTargeting") {
			psi_ = new TargetingMettsType(lrs_,model_,wft_,qn,ioIn_);
		} else if (targeting == "TargetingRixsStatic") {
			psi_ = new TargetingRixsStaticType(lrs_,model_,wft_,qn,ioIn_);
		} else if (targeting == "TargetingRixsDynamic") {
			psi_ = new TargetingRixsDynamicType(lrs_,model_,wft_,qn,ioIn_);
		} else if (targeting == "TargetingExpression") {
			psi_ = new TargetingExpressionType(lrs_, model_, wft_, qn, ioIn_);
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

		VectorStringType targets = {"GroundStateTargeting",
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
		                            "TargetingMultiQ"};

		const SizeType totalTargets = targets.size();
		SizeType count = 0;
		for (SizeType i = 0;i < totalTargets; ++i) {
			if (options.isSet(targets[i])) {
				if (targeting == "AdaptiveDynamicTargeting" &&
				        targets[i] == "DynamicTargeting") continue;
				targeting = targets[i];
				count++;
			}
		}

		if (count == 1) return targeting;

		if (count > 1)
			err("Only one targeting at a time supported\n");

		std::cerr <<" No explicit targeting found, asumming " << targeting <<"\n";
		std::cout <<" No explicit targeting found, asumming " << targeting <<"\n";

		return targeting;
	}

	void check1(PsimagLite::String targeting) const
	{
		if (model_.params().options.isSet("useComplex") &&
		        targeting != "TimeStepTargeting" &&
		        targeting != "ChebyshevTargeting" &&
		        targeting != "GroundStateTargeting" &&
		        targeting != "TargetingCorrelations" &&
		        targeting != "CorrectionTargeting" &&
		        targeting != "CorrectionVectorTargeting" &&
		        targeting != "TargetingInSitu" &&
		        targeting != "TargetingRixsStatic" &&
		        targeting != "TargetingRixsDynamic") {
			err("SolverOptions=useComplex not allowed for " + targeting + "\n");
		}

		if (targeting != "GroundStateTargeting" && BasisType::useSu2Symmetry()) {
			PsimagLite::String str("SU(2) supports only GroundStateTargeting (sorry!)\n");
			throw PsimagLite::RuntimeError(str);
		}
	}

	TargetingBaseType* psi_;
	const LeftRightSuperType& lrs_;
	const ModelType& model_;
	const WaveFunctionTransfType& wft_;
	const typename QnType::VectorQnType& quantumSector_;
	InputValidatorType& ioIn_;
};
}
#endif // TARGETSELECTOR_H

