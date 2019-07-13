#ifndef PREDICATE_AWESOME_H
#define PREDICATE_AWESOME_H

#include "Vector.h"
#include "PredicateAnd.h"

namespace PsimagLite {

/* PSIDOC PredicateAwesome
 Comma-separated list of items, don't leave spaces:
 item0,item1,item2,...

 Each item is either a setting or a predicate

 At least one item is a predicate

 If multiple predicates are present, they are ORed.

 Setting item is of the form
 S1   /\[a-zA-Z]+/
 S2  /\[a-zA-Z]+=\[a-zA-Z\d\.\-]/

 Predicate is a semicolon-separated list of simple predicates,

 SimplePredicate0;SimplePredicate1;...

 where ; means to AND the predicates.

 where SimplePredicate is of the form
 word operator word
 where operator is in {==, <, >, <=, >=, %}
 with the expected meaning, and % means divisible by.
 So l%2 means that the simple predicate is true if l is divisible by 2.

 From this rules, the list of items
 M=2,l%2;l!=0,l=7..11

 will set M=2 and define a predicate that is true if l
 is in the set {2, 4, 6, 7, 8, 9, 10, 12, 14, 16, 18 ...}
 and false otherwise.
 Had we used
 M=2,l%2;l!=0,l=7...11
 the range would have included 11.

 API via SpecType

 S1 calls SpecType::set(String word)
 S2 calls SpecType::set(String word1, String word2)

 The predicate call is
 bool SpecType::isPredicateTrue(const VectorStringType& variables,
								const VectorSizeType& values)
 Only variables of type SizeType are supported for now.

 In the future we could use
 template<typename T>
 bool SpecType::isPredicateTrue(String name, T val)

 template<typename T1, typename T2>
 bool SpecType::isPredicateTrue(String name1, T1 val2, String name2, T2 val2)

 etc.
 */

template<typename SpecType>
class PredicateAwesome {

public:

	typedef Vector<PredicateAnd>::Type VectorPredicateAndType;
	typedef PredicateAnd::VectorStringType VectorStringType;

	PredicateAwesome(String pred)
	    : pred_(pred)
	{
		VectorStringType tokens;
		split(tokens, pred, ",");
		const SizeType n = tokens.size();
		for (SizeType i = 0; i < n; ++i)
			predicateAnd_.push_back(PredicateAnd(tokens[i]));
	}

	template<typename T>
	bool isTrue(String name, T val)
	{
		SizeType n = predicateAnd_.size();
		for (SizeType i = 0; i < n; ++i)
			if (predicateAnd_[i].isTrue(name, val)) return true;
		return false;
	}

private:

	String pred_;
	VectorPredicateAndType predicateAnd_;
};
}
#endif // PREDICATE_AWESOME_H
