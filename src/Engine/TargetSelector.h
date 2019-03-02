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

namespace Dmrg {

template<typename TargetingBaseType>
class TargetSelector {

	typedef typename TargetingBaseType::MatrixVectorType MatrixVectorType;
	typedef typename MatrixVectorType::ModelType ModelType;
	typedef typename ModelType::RealType RealType;
	typedef typename ModelType::ModelHelperType ModelHelperType;
	typedef typename ModelHelperType::LeftRightSuperType LeftRightSuperType;
	typedef typename TargetingBaseType::WaveFunctionTransfType WaveFunctionTransfType;
	typedef typename LeftRightSuperType::BasisWithOperatorsType BasisWithOperatorsType;
	typedef typename BasisWithOperatorsType::BasisType BasisType;
	typedef typename BasisType::QnType QnType;
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

public:

	TargetSelector(const LeftRightSuperType& lrs,
	               const ModelType& model,
	               const WaveFunctionTransfType& wft,
	               const QnType& quantumSector,
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

	TargetingBaseType& operator()(PsimagLite::String targeting)
	{
		if (targeting=="TimeStepTargeting" || targeting == "TargetingAncilla") {
			psi_ = new TargetingTimeStepType(lrs_,model_,wft_,quantumSector_,ioIn_);
		} else if (targeting=="DynamicTargeting") {
			psi_ = new TargetingDynamicType(lrs_,model_,wft_,quantumSector_,ioIn_);
		} else if (targeting=="CorrectionVectorTargeting") {
			psi_ = new TargetingCorrectionVectorType(lrs_,model_,wft_,quantumSector_,ioIn_);
		} else if (targeting=="TargetingChebyshev") {
			psi_ = new TargetingChebyshevType(lrs_,model_,wft_,quantumSector_,ioIn_);
		} else if (targeting=="CorrectionTargeting") {
			psi_ = new TargetingCorrectionType(lrs_,model_,wft_,quantumSector_,ioIn_);
		} else if (targeting == "GroundStateTargeting") {
			psi_ = new TargetingGroundStateType(lrs_,model_,wft_,quantumSector_,ioIn_);
		} else if (targeting == "MettsTargeting") {
			psi_ = new TargetingMettsType(lrs_,model_,wft_,quantumSector_,ioIn_);
		} else if (targeting == "TargetingRixsStatic") {
			psi_ = new TargetingRixsStaticType(lrs_,model_,wft_,quantumSector_,ioIn_);
		} else if (targeting == "TargetingRixsDynamic") {
			psi_ = new TargetingRixsDynamicType(lrs_,model_,wft_,quantumSector_,ioIn_);
		} else {
			err("Unknown targeting " + targeting + "\n");
		}

		psi_->postCtor();

		return *psi_;
	}

private:

	TargetingBaseType* psi_;
	const LeftRightSuperType& lrs_;
	const ModelType& model_;
	const WaveFunctionTransfType& wft_;
	const QnType& quantumSector_;
	InputValidatorType& ioIn_;
};
}
#endif // TARGETSELECTOR_H

