#ifndef AINURSTORE_H
#define AINURSTORE_H
#include "Vector.h"

namespace PsimagLite {

class Store {

public:

	typedef PsimagLite::Vector<String>::Type VectorStringType;

	enum Type {UNKNOWN, STRING, CHAR, SCALAR, VECTOR, MATRIX}; // FUNCTION, GROUP, HASH,

	enum SubType {UNDEFINED, INTEGER, REAL, COMPLEX};

	enum Attribute {NONE, REQUIRED, CONST};

	Store(Type t, SubType s, Attribute a)
	    : type_(t), subType_(s), attr_(a), value_(0)
	{}

	void procDotified(const VectorStringType& dotified,
	                  String right)
	{
		if (dotified.size() == 1) return;
		if (dotified.size() > 2)
			err("More than 2 dots not implemented\n");

		if (dotified[1] == "typeof") {
			setTypeOf(right);
			return;
		}

		if (dotified[1] == "size") {
			// only vectors support size
			if (type_ != Store::VECTOR)
				err(dotified[0] +" is not a vector\n");
			return;
		}

		err("procDotified, unknown " + dotified[1] + "\n");
	}

private:

	void setTypeOf(String s)
	{
		if (s == "string") {
			type_ = STRING;
			return;
		}

		if (s == "char") {
			type_ = CHAR;
			return;
		}

		if (s == "integer" || s == "real" || s == "complex") {
			type_ = SCALAR;
			subType_ = subTypeFromString(s);
			return;
		}

		VectorStringType words;
		split(words, s, ".");
		if (words.size() == 0) return;

		if (words[0] == "vector") {
			type_ = VECTOR;
			if (words.size() > 1)
				subType_ = subTypeFromString(words[1]);
			return;
		}

		if (words[0] == "matrix") {
			type_ = MATRIX;
			if (words.size() > 1)
				subType_ = subTypeFromString(words[1]);
			return;
		}
	}

	static SubType subTypeFromString(String s)
	{
		if (s == "integer") return INTEGER;
		if (s == "real") return REAL;
		if (s == "complex") return COMPLEX;
		return UNDEFINED;
	}

	Type type_;
	SubType subType_;
	Attribute attr_;
	unsigned char *value_;
};
}
#endif // AINURSTORE_H
