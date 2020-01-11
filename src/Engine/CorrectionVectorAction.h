#ifndef CORRECTIONVECTORACTION_H
#define CORRECTIONVECTORACTION_H
#include "Vector.h"
#include "FreqEnum.h"

namespace Dmrg {

template<typename ComplexOrRealType,
         typename TargetParamsType,
         bool isComplex = PsimagLite::IsComplexNumber<ComplexOrRealType>::True>
class CorrectionVectorActionBase {

public:

	typedef typename PsimagLite::Real<ComplexOrRealType>::Type RealType;
	typedef typename PsimagLite::Vector<RealType>::Type VectorRealType;

	enum ActionEnum {ACTION_IMAG, ACTION_REAL};

	CorrectionVectorActionBase(const TargetParamsType& tstStruct,
	                           RealType E0,
	                           const VectorRealType& eigs)
	    : tstStruct_(tstStruct),E0_(E0),eigs_(eigs)
	{}

	ComplexOrRealType operator()(SizeType k) const
	{
		return (tstStruct_.omega().first == PsimagLite::FREQ_REAL) ? actionWhenFreqReal(k)
		                                                           : actionWhenMatsubara(k);
	}

	void setReal() const;

	void setImag() const { action_ = ACTION_IMAG; }

	static bool isValueComplex() { return isComplex; }

protected:

	ComplexOrRealType actionWhenFreqReal(SizeType k) const;

	RealType actionWhenMatsubara(SizeType k) const
	{
		RealType sign = (tstStruct_.type() == 0) ? -1.0 : 1.0;
		RealType wn = tstStruct_.omega().second;
		RealType part1 =  (eigs_[k] - E0_)*sign;
		RealType denom = part1*part1 + wn*wn;
		return (action_ == ACTION_IMAG) ? wn/denom : -part1 / denom;
	}

	const TargetParamsType& tstStruct_;
	RealType E0_;
	const VectorRealType& eigs_;
	mutable ActionEnum action_;
};

template<typename ComplexOrRealType,
         typename TargetParamsType,
         bool isComplex = PsimagLite::IsComplexNumber<ComplexOrRealType>::True>
class CorrectionVectorAction {};

template<typename ComplexOrRealType, typename TargetParamsType>
class CorrectionVectorAction<ComplexOrRealType, TargetParamsType, false> : public
        CorrectionVectorActionBase<ComplexOrRealType, TargetParamsType> {

public:

	typedef CorrectionVectorActionBase<ComplexOrRealType, TargetParamsType> BaseType;
	typedef typename BaseType::RealType RealType;

	CorrectionVectorAction(const TargetParamsType& tstStruct,
	                       RealType E0,
	                       const typename BaseType::VectorRealType& eigs)
	    : BaseType(tstStruct, E0, eigs) {}

	ComplexOrRealType operator()(SizeType k) const
	{
		return (BaseType::tstStruct_.omega().first == PsimagLite::FREQ_REAL)
		        ? actionWhenFreqReal(k) : BaseType::actionWhenMatsubara(k);
	}

	void setReal() const
	{
		BaseType::action_ = BaseType::ACTION_REAL;
	}

private:

	ComplexOrRealType actionWhenFreqReal(SizeType k) const
	{
		RealType sign = (BaseType::tstStruct_.type() == 0) ? -1.0 : 1.0;
		RealType part1 =  (BaseType::eigs_[k] - BaseType::E0_)*sign + BaseType::tstStruct_.omega().second;
		RealType denom = part1*part1 + BaseType::tstStruct_.eta()*BaseType::tstStruct_.eta();
		return (BaseType::action_ == BaseType::ACTION_IMAG) ? BaseType::tstStruct_.eta()/denom :
		                                                      -part1/denom;
	}
};

template<typename ComplexOrRealType, typename TargetParamsType>
class CorrectionVectorAction<ComplexOrRealType, TargetParamsType, true> : public
        CorrectionVectorActionBase<ComplexOrRealType, TargetParamsType>
{

public:

	typedef CorrectionVectorActionBase<ComplexOrRealType, TargetParamsType> BaseType;
	typedef typename BaseType::RealType RealType;

	CorrectionVectorAction(const TargetParamsType& tstStruct,
	                       RealType E0,
	                       const typename BaseType::VectorRealType& eigs)
	    : BaseType(tstStruct, E0, eigs) {}

	ComplexOrRealType operator()(SizeType k) const
	{
		return (BaseType::tstStruct_.omega().first == PsimagLite::FREQ_REAL)
		        ? actionWhenFreqReal(k) : BaseType::actionWhenMatsubara(k);
	}

	void setReal() const
	{
		err("CorrectionVectorSkeleton::Action: Cannot set to real\n");
	}

private:

	ComplexOrRealType actionWhenFreqReal(SizeType k) const
	{
		RealType sign = (BaseType::tstStruct_.type() == 0) ? -1.0 : 1.0;
		RealType part1 =  (BaseType::eigs_[k] - BaseType::E0_)*sign +
		        BaseType::tstStruct_.omega().second;
		return ComplexOrRealType(part1, BaseType::tstStruct_.eta());
	}
};

}
#endif // CORRECTIONVECTORACTION_H
