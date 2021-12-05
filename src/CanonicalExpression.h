#ifndef PSI_CANONICAL_EXPRESSION_H
#define PSI_CANONICAL_EXPRESSION_H
#include "PsimagLite.h"
#include "QuasiCanonical.h"

namespace PsimagLite {

template<typename ItemSpecType>
class AssignAndDestroy {

public:

	typedef typename ItemSpecType::ResultType T;
	typedef typename ItemSpecType::ComplexOrRealType ComplexOrRealType;

	AssignAndDestroy(const T& t) : t_(t) {}

	void assignBackward(T& t2) const
	{
		t2 = t_;
	}

	void assign(const AssignAndDestroy& t2)
	{
		t_ = t2.t_;
	}

	void plusBackward(T& t2) const
	{
		 t2 += t_;
	}

	void multiply(const AssignAndDestroy& t2)
	{
		t_ *= (t2.t_);
	}

	void multiplyScalar(const ComplexOrRealType& scalar)
	{
		t_ *= scalar;
	}

	const bool isEmpty() const { return ItemSpecType::isEmpty(t_); }

	const bool metaEqual(const T& t2) const
	{
		return ItemSpecType::metaEqual(t_, t2);
	}

	static bool metaEqual(const AssignAndDestroy& t1, const AssignAndDestroy& t2)
	{
		return ItemSpecType::metaEqual(t1.t_, t2.t_);
	}

private:

	T t_;
};

template<typename ItemSpecType,
         typename AssignAndDestroyType = AssignAndDestroy<ItemSpecType> >
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

	bool operator()(ResultType& result,
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
			AssignAndDestroyType assignAndDestroy(resultEmpty);// t_ = resultEmpty
			bool isNotEmpty = procCanonicalTerm(assignAndDestroy, quasiCanonical, i, aux);
			if (!isNotEmpty)
				continue;
			if (!assignAndDestroy.metaEqual(result))
				err("CanonicalExpression: metas not equal\n");

			if (i == 0)
				assignAndDestroy.assignBackward(result);  // result = t_;
			else
				assignAndDestroy.plusBackward(result); // result += t_;
		}

		return (ItemSpecType::isEmpty(result)) ? false : true;
	}

private:

	bool procCanonicalTerm(AssignAndDestroyType& assignAndDestroy,
	                       QuasiCanonicalType& quasiCanonical,
	                       SizeType ind,
	                       AuxiliaryType& aux) const
	{
		String termCanonical = quasiCanonical.term(ind);
		ComplexOrRealType factor = 1.0;
		VectorStringType vecStr;
		split(vecStr, termCanonical, "*");

		for (SizeType i = 0; i < vecStr.size(); ++i) {
			procCanonicalFactor(assignAndDestroy, factor, vecStr[i], quasiCanonical, aux);
		}

		if (assignAndDestroy.isEmpty())
			return false;

		assignAndDestroy.multiplyScalar(factor);
		return true;
	}

	void procCanonicalFactor(AssignAndDestroyType& prev,
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

		AssignAndDestroyType term(itemSpec_(termStr, aux));
		if (prev.isEmpty()) {
			prev.assign(term); // prev = term
			return;
		}

		if (!AssignAndDestroyType::metaEqual(term, prev))
			err("CanonicalExpression: metas not equal\n");

		prev.multiply(term); // prev *= term
	}

	const ItemSpecType& itemSpec_;
};
} // namespace PsimagLite
#endif // PSI_CANONICAL_EXPRESSION_H
