#ifndef PSI_CANONICAL_EXPRESSION_H
#define PSI_CANONICAL_EXPRESSION_H
#include "PsimagLite.h"

namespace PsimagLite {

template<typename ItemSpecType>
class CanonicalExpression {

	typedef typename ItemSpecType::ResultType ResultType;
	typedef typename ItemSpecType::ComplexOrRealType ComplexOrRealType;
	typedef typename ItemSpecType::AuxiliaryType AuxiliaryType;
	typedef typename Real<ComplexOrRealType>::Type RealType;
	typedef Vector<String>::Type VectorStringType;

public:

	CanonicalExpression(const ItemSpecType& itemSpec)
	    : itemSpec_(itemSpec)
	{}

	void operator()(ResultType& result,
	                String expr,
	                const ResultType& resultEmpty,
	                AuxiliaryType& aux) const
	{
		// canonical expressions only for now
		// expr --> exprCanonical
		String exprCanonical = expr; // canonicalize here
		VectorStringType vecStr;

		split(vecStr, exprCanonical, "+");

		for (SizeType i = 0; i < vecStr.size(); ++i) {
			ResultType term = resultEmpty;
			procCanonicalTerm(term, vecStr[i], aux);
			if (!ItemSpecType::metaEqual(term, result))
				err("CanonicalExpression: metas not equal\n");

			if (i == 0)
				result = term;
			else
				result += term;
		}

		if (ItemSpecType::isEmpty(result))
			err("CanonicalExpression: expression result is empty\n");
	}

private:

	void procCanonicalTerm(ResultType& term,
	                       String termCanonical,
	                       AuxiliaryType& aux) const
	{
		ComplexOrRealType factor = 1.0;
		VectorStringType vecStr;
		split(vecStr, termCanonical, "*");

		for (SizeType i = 0; i < vecStr.size(); ++i) {
			procCanonicalFactor(term, factor, vecStr[i], aux);
		}

		if (ItemSpecType::isEmpty(term))
			err("CanonicalExpression: term result is empty\n");

		term *= factor;
	}

	void procCanonicalFactor(ResultType& prev,
	                         ComplexOrRealType& factor,
	                         String termStr,
	                         AuxiliaryType& aux) const
	{
		bool isScalar =  isCanonicalScalar(termStr);
		if (isScalar) {
			ComplexOrRealType f = findFactor(termStr,
			                                 false,
			                                 static_cast<ComplexOrRealType*>(0));
			factor *= f;
			return;
		}

		ResultType term = itemSpec_(termStr, aux);
		if (ItemSpecType::isEmpty(prev)) {
			prev = term;
			return;
		}

		if (!ItemSpecType::metaEqual(term, prev))
			err("CanonicalExpression: metas not equal\n");

		prev *= term;
	}

	bool isCanonicalScalar(String termStr) const
	{
		if (termStr.length() == 0)
			err("CanonicalExpression: term must no be empty\n");

		char c = termStr[0];
		bool isDigit = (c >= '0' && c <= '9');
		return (c == '.' || c == '(' || isDigit);
	}

	// Deal with complex here, FIXME
	ComplexOrRealType findFactor(String termStr,
	                             bool prevHadParens,
	                             RealType* dummy) const
	{
		SizeType l = termStr.length();
		if (l == 0)
			err("CanonicalExpression: scalar must no be empty\n");

		char c = termStr[0];
		bool isDigit = (c >= '0' && c <= '9');
		if (c == '.') termStr  = "0" + termStr;
		if (isDigit || c == '.') return atof(termStr.c_str());

		if (c == '-') {
			if (!prevHadParens)
				err("Negative scalar must have enclosing parens\n");

			return atof(termStr.c_str());
		}

		if (c != '(' || termStr[l-1] != ')')
			err("CanonicalExpression: expected enclosing parens\n");

		String tmp = termStr.substr(1, l-2);
		return findFactor(tmp, true, dummy);
	}

	ComplexOrRealType findFactor(String termStr,
	                             bool prevHadParens,
	                             std::complex<RealType>*) const
	{
		std::cerr<<"WARNING: CanonicalExpression: ";
		std::cerr<<"Complex scalars not yet implemented (sorry)\n";
		RealType *x = 0;
		return findFactor(termStr, prevHadParens, x);
	}

	const ItemSpecType& itemSpec_;
};
} // namespace PsimagLite
#endif // PSI_CANONICAL_EXPRESSION_H
