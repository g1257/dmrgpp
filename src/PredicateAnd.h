#ifndef PREDICATEAND_H
#define PREDICATEAND_H
#include "Vector.h"
#include "PredicateSimple.h"
#include "PsimagLite.h"

namespace Dmrg {

/* PSIDOC PredicateAnd
 PredicateAnd is a semicolon-separated list of simple predicates,

 SimplePredicate0;SimplePredicate1;...

 where ; means to AND the predicates.
 */
class PredicateAnd {

public:

	typedef PsimagLite::Vector<PredicateSimple>::Type VectorPredicateSimpleType;

	PredicateAnd(PsimagLite::String pred)
	    : pred_(pred)
	{
		PsimagLite::split(tokens, str, ";");
		for (SizeType i = 0; i < n; ++i)
			vPredicateSimple_.push_back(PredicateSimple(str[i]));
	}

	template<typename T>
    bool isTrue(PsimagLite::String name, T val)
	{
		SizeType n = vPredicateSimple_.size();
		for (SizeType i = 0; i < n; ++i)
			if (!vPredicateSimple_[i].isTrue(name, val)) return false;
		return true;
	}

private:

	PsimagLite::String pred_;
	VectorPredicateSimpleType vPredicateSimple_;
};
}
#endif // PREDICATEAND_H
