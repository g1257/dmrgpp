#ifndef PSI_CANONICAL_EXPRESSION_H
#define PSI_CANONICAL_EXPRESSION_H
#include "PsimagLite.h"
#include "QuasiCanonical.h"

namespace PsimagLite {

template<typename ItemSpecType>
class CanonicalExpression {

	typedef typename ItemSpecType::ResultType ResultType;
	typedef typename ItemSpecType::ComplexOrRealType ComplexOrRealType;
	typedef typename ItemSpecType::AuxiliaryType AuxiliaryType;
	typedef typename Real<ComplexOrRealType>::Type RealType;
	typedef Vector<String>::Type VectorStringType;
	typedef PsimagLite::QuasiCanonical<ComplexOrRealType> QuasiCanonicalType;

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
		QuasiCanonicalType quasiCanonical(expr);
		//String exprCanonical = expr; // canonicalize here
		//VectorStringType vecStr;
		const SizeType numberOfTerms = quasiCanonical.numberOfTerms();

		for (SizeType i = 0; i < numberOfTerms; ++i) {
			ResultType term = resultEmpty;
			procCanonicalTerm(term, quasiCanonical, i, aux);
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
	                       QuasiCanonicalType& quasiCanonical,
	                       SizeType ind,
	                       AuxiliaryType& aux) const
	{
		String termCanonical = quasiCanonical.term(ind);
		ComplexOrRealType factor = 1.0;
		VectorStringType vecStr;
		split(vecStr, termCanonical, "*");

		for (SizeType i = 0; i < vecStr.size(); ++i) {
			procCanonicalFactor(term, factor, vecStr[i], quasiCanonical, aux);
		}

		if (ItemSpecType::isEmpty(term))
			err("CanonicalExpression: term result is empty\n");

		term *= factor;
	}

	void procCanonicalFactor(ResultType& prev,
	                         ComplexOrRealType& factor,
	                         String termStr,
	                         QuasiCanonicalType& quasiCanonical,
	                         AuxiliaryType& aux) const
	{
		bool isCaScalar = QuasiCanonicalType::isRealScalar(termStr);
		if (isCaScalar) {
			const ComplexOrRealType f = PsimagLite::atof(termStr);
			factor *= f;
			return;
		}

		bool isPureCmplx = QuasiCanonicalType::isPureComplex(termStr);
		if (isPureCmplx) {
			const ComplexOrRealType f = QuasiCanonicalType::pureComplex(termStr);
			factor *= f;
			return;
		}

		const int ind = quasiCanonical.scalarIndex(termStr);
		if (ind >= 0) {
			const ComplexOrRealType f = quasiCanonical.scalarFromIndex(ind);
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

	const ItemSpecType& itemSpec_;
};
} // namespace PsimagLite
#endif // PSI_CANONICAL_EXPRESSION_H
