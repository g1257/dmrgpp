#ifndef PREDICATE_SIMPLE_H
#define PREDICATE_SIMPLE_H
#include "Vector.h"
#include "PsimagLite.h"
#include "ExpressionCalculator.h"

/* PSIDOC PredicateSimple
 PredicateSimple is of the form
 word operator word
 where operator is in {==, <, >, <=, >=, %%}
 with the usual meaning, and %% means divisible by.
 So l%%2 means that the simple predicate is true if l is divisible by 2.
 */
namespace PsimagLite {

class PredicateSimple {

public:

	typedef Vector<String>::Type VectorStringType;

	PredicateSimple(String pred)
	    : pred_(pred)
	{
		SizeType length = 0;
		size_t location = String::npos;

		SizeType n = pred.length();
		if (n < 3) err("PredicateSimple: pred must have at least 3 characters\n");

		for (SizeType i = 0; i < n - 1; ++i) {

			// order matters: test with LONGER FIRST
			String maybeOp = pred.substr(i, 2);
			if (operatorLength(maybeOp) == 2) {
				location = i;
				length = 2;
				break;
			}

			maybeOp = pred.substr(i, 1);
			if (operatorLength(maybeOp) == 1) {
				location = i;
				length = 1;
				break;
			}
		}

		if (length == 0)
			err("Could not find operator in predicate " + pred + "\n");

		lhs_ = pred.substr(0, location);
		op_ = pred.substr(location, length);
		rhs_ = pred.substr(location + length, n - location - length);
		if (lhs_ == "" || rhs_ == "")
			err("Left or right expression is empty\n");
	}

	template<typename SomeVectorType>
	bool isTrue(const VectorStringType& names, const SomeVectorType& vals)
	{
		typedef typename SomeVectorType::value_type SomeValueType;
		SomeValueType lv = getValue(lhs_, names, vals);
		SomeValueType rv = getValue(rhs_, names, vals);
		return compareOnOp(lv, op_, rv);
	}

private:

	template<typename T>
	static bool compareOnOp(T lv, String op, T rv)
	{
		// {==, <, >, <=, >=, %%}
		// If you add something here, add it also to PredicateSimple.cpp
		if (op == "==")
			return (lv == rv);
		if (op == "!=")
			return (lv != rv);
		if (op == "<")
			return (lv < rv);
		if (op == ">")
			return (lv > rv);
		if (op == "<=")
			return (lv <= rv);
		if (op == ">=")
			return (lv >= rv);
		if (op == "%%")
			return ((lv % rv) == 0);
		throw RuntimeError("Unknown operator " + op + "\n");
	}

	static SizeType operatorLength(String op)
	{
		const SizeType n = ops_.size();
		for (SizeType i = 0; i < n; ++i)
			if (op == ops_[i]) return op.length();

		return 0;
	}

	template<typename SomeVectorType>
	static typename SomeVectorType::value_type getValue(String hs,
	                                                    const VectorStringType& names,
	                                                    const SomeVectorType& vals)
	{
		String numericHs = replaceVariables(hs, names, vals);
		VectorStringType tokens;
		split(tokens, numericHs, "|");
		typedef ExpressionCalculator<typename SomeVectorType::value_type> ExpressionCalculatorType;
		ExpressionCalculatorType expressionCalculator(tokens);
		return expressionCalculator();
	}

	template<typename SomeVectorType>
	static String replaceVariables(String hs,
	                               const VectorStringType& names,
	                               const SomeVectorType& vals)
	{
		const SizeType n = names.size();
		assert(n == vals.size());
		String buffer = hs;
		for (SizeType i = 0; i < n; ++i) {
			buffer = replaceVariable(buffer, names[i], vals[i]);
		}

		return buffer;
	}

	template<typename T>
	static String replaceVariable(String hs,
	                              String name,
	                              T val)
	{
		const String valString = ttos(val);
		const SizeType nameLength = name.length();

		while (true) {
			size_t index = hs.find(name);
			if (index == String::npos) return hs;
			String part1 = (index == 0) ? "" : hs.substr(0, index);
			String part2 = hs.substr(index + nameLength, hs.length() - nameLength - index);
			hs = part1 + valString + part2;
		}

		return hs;
	}

	static VectorStringType ops_;
	String pred_;
	String lhs_;
	String op_;
	String rhs_;
};
}
#endif // PREDICATE_SIMPLE_H
