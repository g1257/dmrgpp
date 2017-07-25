#ifndef AINURSTORE_H
#define AINURSTORE_H
#include "Vector.h"

namespace PsimagLite {

class Store {

public:

	typedef PsimagLite::Vector<String>::Type VectorStringType;

	enum Type {UNKNOWN, STRING, CHAR, SCALAR, VECTOR, MATRIX, FUNCTION}; // GROUP, HASH,

	enum SubType {UNDEFINED, INTEGER, REAL, COMPLEX};

	enum Attribute {NONE, REQUIRED, CONST};

//	Store(Type t, SubType s, Attribute a)
//	    : type_(t), subType_(s), attr_(a), value_(0)
//	{}

	Store(String s)
	    : type_(UNKNOWN), subType_(UNDEFINED), attr_(NONE), used_(0)
	{
		setTypeOf(s);
	}

	Store(String s, String a)
	    : type_(UNKNOWN), subType_(UNDEFINED), attr_(NONE), used_(0)
	{
		setTypeOf(s);
		setAttr(a);
	}

	void setRhs(String rhs)
	{
		value_.clear();
		switch (type_) {
		case STRING:
		case CHAR:
		case FUNCTION:
		case SCALAR:
			value_.push_back(rhs);
			break;
		case VECTOR:
			setVectorValue(rhs);
			break;
		case MATRIX:
			setMatrixValue(rhs);
			break;
		default:
			break;
		}
	}

	Type type() const { return type_; }

	SubType subType() const { return subType_; }

	String value(SizeType ind) const
	{
		assert(ind < value_.size());
		return value_[ind];
	}

	void increaseUsage() const { ++used_; }

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

	void setVectorValue(String rhs)
	{
		SizeType last = rhs.size();
		if (last == 0) return;
		--last;
		if (rhs[0] != '[' || rhs[last] != ']')
			err("Malformed vector " + rhs + "\n");
		split(value_, rhs, ",");
	}

	void setMatrixValue(String rhs)
	{
		std::cerr<<"Not setting matrix value to "<<rhs<<"\n";
	}

	Type type_;
	SubType subType_;
	Attribute attr_;
	VectorStringType value_;
	mutable SizeType used_;
};
}
#endif // AINURSTORE_H
