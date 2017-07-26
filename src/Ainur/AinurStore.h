#ifndef AINURSTORE_H
#define AINURSTORE_H
#include "Vector.h"

namespace PsimagLite {

class Store {

public:

	typedef PsimagLite::Vector<String>::Type VectorStringType;

	enum Type {UNKNOWN, SCALAR, VECTOR, MATRIX}; // HASH, FUNCTION

	enum SubType {UNDEFINED, INTEGER, REAL, COMPLEX, STRING, CHAR, GROUP};

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
			std::cerr<<"setRhs not implemented, rhs= "<<rhs<<"\n";
			break;
		}
	}

	Type type() const { return type_; }

	SubType subType() const { return subType_; }

	SizeType valueSize() const { return value_.size(); }

	String value(SizeType ind) const
	{
		if (ind >= value_.size())
			throw RuntimeError("Not such value\n");

		return value_[ind];
	}

	void increaseUsage() const { ++used_; }

	SizeType used() const { return used_; }

private:

	void setTypeOf(String s)
	{
		subType_ = subTypeFromString(s);
		if (subType_ != UNDEFINED) {
			type_ = SCALAR;
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
		if (s == "char") return CHAR;
		if (s == "string") return STRING;
		if (s == "group") return GROUP;
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
		SizeType last = rhs.length();
		if (last < 2)
			err("Vector must be enclosed in brakets\n");

		--last;
		if (rhs[0] != '[' || rhs[last] != ']')
			err("Malformed vector " + rhs + "\n");
		rhs = (last < 2 ) ? "" : rhs.substr(1,last - 1);
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
