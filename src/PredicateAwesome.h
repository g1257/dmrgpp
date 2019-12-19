#ifndef PREDICATE_AWESOME_H
#define PREDICATE_AWESOME_H

#include "Vector.h"
#include "PredicateAnd.h"
#include "PredicateDefaultSpec.h"

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
 where operator is in {==, <, >, <=, >=, %%}
 with the expected meaning, and %% means divisible by.
 So l%%2 means that the simple predicate is true if l is divisible by 2.

 From this rules, the list of items
 M=2,l%%2;l!=0,l>=7,l<11
 will set M=2 and define a predicate that is true if l
 is in the set {2, 4, 6, 7, 8, 9, 10, 12, 14, 16, 18 ...}
 and false otherwise.

 API via SpecType

 S1 calls SpecType::set(String word)
 S2 calls SpecType::set(String word1, String word2)

 The predicate call is

 template<typename T>
 bool isPredicateTrue(String name, T val)

 template<typename T1, typename T2>
 bool isPredicateTrue(String name1, T1 val2, String name2, T2 val2)

 etc. Can we do this with varargs templates?? FIXME TODO
 */

template<typename SpecType = PredicateDefaultSpec>
class PredicateAwesome {

public:

	typedef Vector<PredicateAnd>::Type VectorPredicateAndType;
	typedef PredicateAnd::VectorStringType VectorStringType;

	PredicateAwesome(String pred, char orSep = ',')
	    : pred_(pred)
	{
		if (pred_ == "") return;
		VectorStringType tokens;
		String orSeparator(",");
		orSeparator[0] = orSep;
		split(tokens, pred, orSeparator);
		const SizeType n = tokens.size();
		for (SizeType i = 0; i < n; ++i)
			predicateAnd_.push_back(PredicateAnd(tokens[i]));
	}

	template<typename T>
	bool isTrue(String name, T val)
	{
		if (pred_ == "") return false;
		SizeType n = predicateAnd_.size();
		for (SizeType i = 0; i < n; ++i)
			if (predicateAnd_[i].isTrue(name, val)) return true;
		return false;
	}

	template<typename T1, typename T2>
	bool isTrue(String name1, T1 val1, String name2, T2 val2)
	{
		if (pred_ == "") return false;
		SizeType n = predicateAnd_.size();
		for (SizeType i = 0; i < n; ++i)
			if (predicateAnd_[i].isTrue(name1, val1, name2, val2)) return true;
		return false;
	}

	template<typename T1, typename T2>
	bool isTrue(String name1, T1 val1,
	            String name2, T2 val2,
	            String name3, T1 val3)
	{
		if (pred_ == "") return false;
		SizeType n = predicateAnd_.size();
		for (SizeType i = 0; i < n; ++i)
			if (predicateAnd_[i].isTrue(name1, val1,
			                            name2, val2,
			                            name3, val3)) return true;
		return false;
	}

	template<typename T1, typename T2>
	bool isTrue(String name1, T1 val1,
	            String name2, T2 val2,
	            String name3, T1 val3,
	            String name4, T2 val4)
	{
		if (pred_ == "") return false;
		SizeType n = predicateAnd_.size();
		for (SizeType i = 0; i < n; ++i)
			if (predicateAnd_[i].isTrue(name1, val1,
			                            name2, val2,
			                            name3, val3,
			                            name4, val4)) return true;
		return false;
	}

	static void replaceAll(String& str, const String& from, const String& to)
	{
		if (from.empty()) return;

		size_t start_pos = 0;

		while ((start_pos = str.find(from, start_pos)) != String::npos) {
			str.replace(start_pos, from.length(), to);
			start_pos += to.length(); // In case 'to' contains 'from', like replacing 'x' with 'yx'
		}
	}

private:

	String pred_;
	VectorPredicateAndType predicateAnd_;
};
}
#endif // PREDICATE_AWESOME_H
