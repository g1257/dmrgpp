#ifndef RAHULOPERATOR_H
#define RAHULOPERATOR_H
#include "Vector.h"

namespace LanczosPlusPlus {

template <typename ComplexOrRealType> class RahulOperator {

	enum class Label
	{
		IDENTITY,
		N,
		SZ,
		C
	};

public:

	RahulOperator(PsimagLite::String label, SizeType dof, bool transpose = false)
	    : dof_(dof)
	    , transpose_(transpose)
	{
		label_ = fromString(label);
	}

	SizeType dof() const { return dof_; }

	// The output of this function should be interpreted as follows.
	// If the function returns false, bit and result MUST BE ignored
	// If the function results true, bit is what the bit should be and result
	// contains the multiplier
	bool actOn(bool& bit, ComplexOrRealType& result) const
	{
		static const ComplexOrRealType zeroPointFive = 0.5;
		result                                       = 1;
		const bool bitSaved                          = bit;
		switch (label_) {
		case Label::IDENTITY:
			return true;
			break;
		case Label::N:
			return (bitSaved);
			break;
		case Label::SZ:
			result = (bitSaved) ? -zeroPointFive : zeroPointFive;
			return true;
			break;
		case Label::C:
			bit = !bit;
			return ((bitSaved && !transpose_) || (!bitSaved && transpose_));
			break;
		default:
			throw PsimagLite::RuntimeError("RahulOperator::actOn internal error\n");
			break;
		}
	}

	bool isFermionic() const { return (label_ == Label::C); }

private:

	static Label fromString(PsimagLite::String l)
	{
		if (l == "c")
			return Label::C;
		if (l == "identity")
			return Label::IDENTITY;
		if (l == "sz")
			return Label::SZ;
		if (l == "n")
			return Label::N;
		throw PsimagLite::RuntimeError("RahulOperator: Unknow label " + l + "\n");
	}

	SizeType dof_;
	bool     transpose_;
	Label    label_;
};
}
#endif // RAHULOPERATOR_H
