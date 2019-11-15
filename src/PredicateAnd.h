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
		VectorStringType names{name};
		typename Vector<T>::Type values{val};
		SizeType n = vPredicateSimple_.size();
		for (SizeType i = 0; i < n; ++i)
			if (!vPredicateSimple_[i].isTrue(names, values)) return false;
		return true;
	}

	template<typename T>
	bool isTrue(String name1, T val1, String name2, T val2)
	{
		VectorStringType names{name1, name2};
		typename Vector<T>::Type values{val1, val2};
		SizeType n = vPredicateSimple_.size();
		for (SizeType i = 0; i < n; ++i)
			if (!vPredicateSimple_[i].isTrue(names, values)) return false;
		return true;
	}

	template<typename T>
	bool isTrue(String name1, T val1,
	            String name2, T val2,
	            String name3, T val3,
	            String name4, T val4)
	{
		VectorStringType names{name1, name2, name3, name4};
		typename Vector<T>::Type values{val1, val2, val3, val4};
		SizeType n = vPredicateSimple_.size();
		for (SizeType i = 0; i < n; ++i)
			if (!vPredicateSimple_[i].isTrue(names, values)) return false;
		return true;
	}

private:

	String pred_;
	VectorPredicateSimpleType vPredicateSimple_;
};
}
#endif // PREDICATEAND_H
