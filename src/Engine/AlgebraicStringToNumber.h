#ifndef ALGEBRAICSTRINGTONUMBER_H
#define ALGEBRAICSTRINGTONUMBER_H
#include "CanonicalExpression.h"
#include "Vector.h"

namespace Dmrg {

template <typename FieldType> class AlgebraicStringToNumber {

	class LoopLengthSpec {

	public:

		typedef FieldType ResultType;
		typedef FieldType ComplexOrRealType;
		typedef FieldType AuxiliaryType;

		LoopLengthSpec(SizeType n)
		    : n_(n)
		    , c_(n / 2)
		{
			if (n % 2 == 1)
				c_ = 0;
		}

		static bool isEmpty(ResultType x) { return (x == 0); }

		static bool metaEqual(ResultType, ResultType) { return true; }

		ResultType operator()(PsimagLite::String str, FieldType numberOfSites) const
		{
			if (str == "%lh")
				return n_ - 2;

			if (str == "%n")
				return n_;

			if (str == "%c") {
				if (c_ == 0)
					err("%c cannot be used with odd number of sites\n");
				return c_;
			}

			return PsimagLite::atof(str);
		}

		PsimagLite::String convertVal(PsimagLite::String val) const
		{
			if (val == "-%lh")
				return "(-1.0)*%lh";
			return val;
		}

	private:

		SizeType n_;
		SizeType c_;
	};

public:

	AlgebraicStringToNumber(PsimagLite::String msg, SizeType numberOfSites)
	    : msg_(msg)
	    , numberOfSites_(numberOfSites)
	{ }

	int procLength(PsimagLite::String val)
	{
		LoopLengthSpec loopLengthSpec(numberOfSites_);

		val = loopLengthSpec.convertVal(val);

		FieldType numberOfSites = numberOfSites_;
		PsimagLite::CanonicalExpression<LoopLengthSpec> canonicalExpression(loopLengthSpec);
		typename LoopLengthSpec::ResultType opEmpty(1);
		typename LoopLengthSpec::ResultType p(1);
		canonicalExpression(p, val, opEmpty, numberOfSites);
		if (p != static_cast<int>(p))
			err(msg_ + ": string value " + val + " yields non integer\n");
		return p;
	}

private:

	PsimagLite::String msg_;
	SizeType numberOfSites_;
};
}
#endif // ALGEBRAICSTRINGTONUMBER_H
