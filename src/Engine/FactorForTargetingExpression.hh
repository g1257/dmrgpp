#ifndef FACTORFORTARGETINGEXPRESSION_HH
#define FACTORFORTARGETINGEXPRESSION_HH
#include "Vector.h"
#include <string>

namespace Dmrg {

template<typename ComplexOrRealType>
class FactorForTargetingExpression {

	using RealType = typename PsimagLite::Real<ComplexOrRealType>::Type;
	using VectorType = std::vector<ComplexOrRealType>;
	using VecStringType = std::vector<std::string>;

public:
	FactorForTargetingExpression(const ComplexOrRealType& value)
	{
		factors_.push_back(value);
		strFactors_.push_back("");
		setStr(0);
	}

	void multiply(const ComplexOrRealType& val)
	{
		if (factors_.size() != 1) {
			err("FactorForTargetingExpression: Cannot multiply vector by scalar\n");
		}


                factors_[0] *= val;
		setStr(0);
	}

	void set(const ComplexOrRealType& val)
	{
		factors_.resize(1);

		factors_[0] = val;
		setStr();
	}

	void push(const ComplexOrRealType& val)
	{
		SizeType index = factors_.size();
		factors_.push_back(val);
		strFactors_.push_back("");
		setStr(index);
	}

	std::string toString() const
	{
		SizeType n = strFactors_.size();
		if (n == 0) {
			return "";
		}

		if (n == 1) {
			return (strFactors_[0] != "") ? strFactors_[0] + "*" : "";
		}

		std::string buffer("[");
		buffer += strFactors_[0];
		for (SizeType i = 1; i < n; ++i) {
			buffer += "," + strFactors_[i];
		}

		buffer += "]*";
		return buffer;
	}

	ComplexOrRealType value() const
	{
		if (factors_.size() != 1) {
			err("FactorForTargetingExpression: Cannot convert vector to string\n");
		}

		return factors_[0];
	}

	ComplexOrRealType value(SizeType ind) const
	{
		assert(ind < factors_.size());
		return factors_[ind];
	}

	void swap()
	{
		if (factors_.size() != 2)
			err("Cannot call swap without 2 factors\n");

		ComplexOrRealType f = factors_[0];
		factors_[0] = factors_[1];
		factors_[1] = f;
		setStr();
	}

private:

	void setStr()
	{
		strFactors_.resize(factors_.size());
		unsigned int n = strFactors_.size();
		for (SizeType i = 0; i < n; ++i) {
			setStr(i);
		}
	}

	void setStr(SizeType ind)
	{
		return (PsimagLite::IsComplexNumber<ComplexOrRealType>::True) ? setStrImagOne(ind)
									      : setStrRealOne(ind);
	}

	void setStrRealOne(SizeType ind)
	{
		assert(factors_.size() == strFactors_.size());
		assert(ind < factors_.size());

		assert(PsimagLite::imag(factors_[ind]) == 0);

		const RealType f = PsimagLite::real(factors_[ind]);

		if (f < 0)
			strFactors_[ind] = "(" + ttos(f) + ")";
		else
			strFactors_[ind] = ttos(f);

		if (f == 1)
			strFactors_[ind] = "";
	}

	void setStrImagOne(SizeType ind)
	{
		const RealType freal = PsimagLite::real(factors_[ind]);
		const RealType fimag = PsimagLite::imag(factors_[ind]);

		strFactors_[ind] = "(" + ttos(freal) + "+" + ttos(fimag) + "i)";

		if (freal == 1 && fimag == 0)
			strFactors_[ind] = "";
	}

	VectorType factors_;
	VecStringType strFactors_;
};
}
#endif // FACTORFORTARGETINGEXPRESSION_HH
