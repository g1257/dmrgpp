#ifndef SPECFORTARGETINGEXPRESSION_H
#define SPECFORTARGETINGEXPRESSION_H
#include "Vector.h"

namespace Dmrg {

template<typename VectorWithOffsetType>
class AlgebraForTargetingExpression {

public:

	typedef typename VectorWithOffsetType::value_type ComplexOrRealType;
	typedef int AuxiliaryType;

	const VectorWithOffsetType& toVectorWithOffset() const { return data_; }

	AlgebraForTargetingExpression() : finalized_(false), factor_(1.0) {}

	AlgebraForTargetingExpression(PsimagLite::String str)
	    : finalized_(false), str_(str), factor_(1.0) {}

	// AlgebraForTargetingExpression(const AlgebraForTargetingExpression&) = delete;

	// AlgebraForTargetingExpression& operator=(const AlgebraForTargetingExpression&) = delete;

	AlgebraForTargetingExpression& operator+=(const AlgebraForTargetingExpression& other)
	{
		AlgebraForTargetingExpression otherCopy = other;
		otherCopy.finalize();
		finalize();
		data_ += other.toVectorWithOffset();
		return *this;
	}

	AlgebraForTargetingExpression& operator*=(const AlgebraForTargetingExpression& other)
	{
		if (other.finalized_ || finalized_)
			err("AlgebraForTargetingExpression: Two finalized terms cannot be multiplied\n");
		if (other.str_ == "") return *this;

		str_ += "*";
		str_ += other.str_;
		return *this;
	}

	AlgebraForTargetingExpression& operator*=(const ComplexOrRealType& scalar)
	{
		factor_ *= scalar;
		return *this;
	}

	const PsimagLite::String& toString() const { return str_; }

	bool finalized() const { return finalized_; }

private:

	void finalize()
	{
		if (finalized_) return;
		finalized_ = true;
		// one str_ must be |something>
		// get |something>
		// loop over str_ other than |something> and apply
		err("AlgebraForTargetingExpression::finalize() not implemented\n");
	}

	bool finalized_;
	PsimagLite::String str_;
	VectorWithOffsetType data_;
	ComplexOrRealType factor_;
};

template<typename VectorWithOffsetType>
class SpecForTargetingExpression {

public:

	typedef AlgebraForTargetingExpression<VectorWithOffsetType> AlgebraType;
	typedef AlgebraType ResultType;
	typedef typename VectorWithOffsetType::value_type ComplexOrRealType;
	typedef typename AlgebraType::AuxiliaryType AuxiliaryType;

	AlgebraType operator()(PsimagLite::String str, AuxiliaryType) const
	{
		return AlgebraType(str);
	}

	static bool metaEqual(const ResultType&, const ResultType&) { return true; }

	static bool isEmpty(const ResultType& term)
	{
		if (term.finalized()) return false;
		return (term.toString() == "");
	}
};
}
#endif // SPECFORTARGETINGEXPRESSION_H
