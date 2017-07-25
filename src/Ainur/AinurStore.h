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

//	Store(Type t, SubType s, Attribute a)
//	    : type_(t), subType_(s), attr_(a), value_(0)
//	{}

	Store(String s): type_(UNKNOWN), subType_(UNDEFINED), attr_(NONE)
	{
		setTypeOf(s);
	}

	Store(String s, String a) : type_(UNKNOWN), subType_(UNDEFINED), attr_(NONE)
	{
		setTypeOf(s);
		setAttr(a);
	}

	void setRhs(String rhs)
	{
		std::cerr<<"Not setting RHS to "<<rhs<<"\n";
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

	void setAttr(String s)
	{
		if (s == "require") {
			attr_ = REQUIRED;
			return;
		}

		if (s == "const") {
			attr_ = CONST;
			return;
		}

		err("Unknown attribute " + s + "\n");
	}

	Type type_;
	SubType subType_;
	Attribute attr_;
};
}
#endif // AINURSTORE_H
