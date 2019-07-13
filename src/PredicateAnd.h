#ifndef PREDICATEAND_H
#define PREDICATEAND_H
#include "Vector.h"
#include "PredicateSimple.h"
#include "PsimagLite.h"

namespace PsimagLite {

/* PSIDOC PredicateAnd
 PredicateAnd is a semicolon-separated list of simple predicates,

 SimplePredicate0;SimplePredicate1;...

 where ; means to AND the predicates.
 */
class PredicateAnd {

public:

	typedef Vector<PredicateSimple>::Type VectorPredicateSimpleType;
	typedef PredicateSimple::VectorStringType VectorStringType;

	PredicateAnd(String pred)
	    : pred_(pred)
	{
		VectorStringType tokens;
		split(tokens, pred, ";");
		const SizeType n = tokens.size();
		for (SizeType i = 0; i < n; ++i)
			vPredicateSimple_.push_back(PredicateSimple(tokens[i]));
	}

	template<typename T>
	bool isTrue(String name, T val)
	{
		SizeType n = vPredicateSimple_.size();
		for (SizeType i = 0; i < n; ++i)
			if (!vPredicateSimple_[i].isTrue(name, val)) return false;
		return true;
	}

private:

	String pred_;
	VectorPredicateSimpleType vPredicateSimple_;
};
}
#endif // PREDICATEAND_H
