#ifndef QUASICANONICAL_H
#define QUASICANONICAL_H
#include "AST/ExpressionForAST.h"
#include "AST/PlusMinusMultiplyDivide.h"
#include "CmplxOrReal.h"
#include "PsimagLite.h"
#include "Vector.h"

namespace PsimagLite
{

// level one parens only
// TODO FIXME : Generalize to multiple levels
template <typename ComplexOrRealType>
class QuasiCanonical
{

public:

	typedef Vector<String>::Type VectorStringType;
	typedef Vector<bool>::Type VectorBoolType;
	typedef typename Vector<ComplexOrRealType>::Type VectorType;
	typedef typename Real<ComplexOrRealType>::Type RealType;
	typedef std::complex<RealType> ComplexType;

	QuasiCanonical(String str)
	    : str_(str)
	{
		const SizeType len = str_.length();
		String tempBuffer;
		String status = "closed";
		char prev = '\0';
		for (SizeType i = 0; i < len; ++i) {

			if (str_[i] == '-' && i > 0) {
				if (prev != '(' && prev != '+') {
					throw RuntimeError(
					    "The - sign must be preceded by "
					    "nothing, parens, or +\n");
				}
			}

			prev = str_[i];

			if (str_[i] == '(') {
				if (status == "open")
					throw RuntimeError(
					    "Too many parens levels (one only "
					    "supported for now)\n");
				status = "open";
				continue;
			}

			if (str_[i] == ')') {
				if (status == "closed")
					throw RuntimeError(
					    "Unbalanced parens, closed\n");
				status = "closed";
				mainBuffer_ += "@" + ttos(ats_.size()) + "@";
				ats_.push_back(tempBuffer);
				tempBuffer = "";
				continue;
			}

			if (status == "closed") {
				mainBuffer_ += str_[i];
			} else {
				tempBuffer += str_[i];
			}
		}

		if (status == "open")
			throw RuntimeError("Unbalanced parens, open\n");

		split(terms_, mainBuffer_, "+");

		const SizeType nscalars = ats_.size();
		cachedValues_.resize(nscalars);
		for (SizeType i = 0; i < nscalars; ++i) {
			cachedValues_[i] = simpleArithmetics(ats_[i]);
		}
	}

	SizeType numberOfTerms() const { return terms_.size(); }

	const String& term(SizeType i) const
	{
		assert(i < terms_.size());
		return terms_[i];
	}

	int scalarIndex(String str) const
	{
		const SizeType len = str.length();
		if (len < 3)
			return -1;
		if (str[0] != '@' || str[len - 1] != '@')
			return -1;

		String snumber = str.substr(1, len - 2);
		SizeType number = atoi(snumber);
		if (number >= ats_.size())
			return -1;
		return number;
	}

	ComplexOrRealType scalarFromIndex(SizeType ind) const
	{
		assert(ind < cachedValues_.size());

		return cachedValues_[ind];
	}

	static bool isPureComplex(String t)
	{
		if (t == "i")
			return true;

		const SizeType n = t.length();
		if (n < 2)
			return false;
		String tmp = t.substr(0, n - 1);
		return isRealScalar(tmp);
	}

	static ComplexOrRealType pureComplex(String t)
	{
		static const bool isComplex = IsComplexNumber<ComplexOrRealType>::True;
		if (!isComplex)
			err("i = sqrt(-1) found in code path that is real\n");

		CpmlxOrReal<RealType, (isComplex) ? 1 : 0> cmplxOrReal(t);
		return cmplxOrReal.value();
	}

	static bool isRealScalar(String termStr)
	{
		const SizeType n = termStr.length();
		if (n == 0)
			err("CanonicalExpression: term must not be empty\n");

		for (SizeType i = 0; i < n; ++i) {
			char c = termStr[i];
			bool isDigit = (c >= '0' && c <= '9');
			if (c == '.' || c == '-' || c == '+' || isDigit)
				continue;
			return false;
		}

		return true;
	}

private:

	ComplexOrRealType simpleArithmetics(String str)
	{
		replaceAll(str, "pi", ttos(M_PI));

		VectorStringType ve;
		split(ve, str, ":");

		typedef PlusMinusMultiplyDivide<ComplexOrRealType>
		    PrimitivesType;
		PrimitivesType primitives;
		ExpressionForAST<PrimitivesType> expresionForAST(ve,
		    primitives);
		return expresionForAST.exec();
	}

	String str_;
	String mainBuffer_;
	VectorStringType ats_;
	VectorStringType terms_;
	VectorType cachedValues_;
};
} // namespace PsimagLite
#endif // QUASICANONICAL_H
