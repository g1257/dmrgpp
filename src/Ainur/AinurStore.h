#ifndef AINURSTORE_H
#define AINURSTORE_H
#include "Vector.h"
#include "AinurLexical.h"

namespace PsimagLite {

class Store {

public:

	typedef AinurLexical AinurLexicalType;
	typedef AinurLexicalType::VectorStringType VectorStringType;

	enum Type {UNKNOWN, SCALAR, VECTOR, MATRIX}; // HASH, FUNCTION

	enum SubType {UNDEFINED, INTEGER, REAL, COMPLEX, STRING, CHAR, GROUP};

	enum Attribute {NONE, REQUIRED, CONST};

	Store(String s, String a)
	    : type_(UNKNOWN), subType_(UNDEFINED), attr_(NONE), used_(0)
	{
		setTypeOf(s);
		if (a != "")
			setAttr(a);
	}

	void setRhs(String rhs, String name)
	{
		value_.clear();
		switch (type_) {
		case SCALAR:
			value_.push_back(rhs);
			break;
		case VECTOR:
			setVectorValue(value_, rhs, name);
			break;
		case MATRIX:
			setMatrixValue(rhs, name);
			break;
		default:
			std::cerr<<"setRhs not implemented, rhs= "<<rhs<<"\n";
			break;
		}
	}

	Type type() const { return type_; }

	SubType subType() const { return subType_; }

	SizeType valueSize() const { return value_.size(); }

	const String& value(SizeType ind, String name) const
	{
		if (ind >= value_.size())
			throw RuntimeError("No value for " + name + "\n");

		return value_[ind];
	}

	String& value(SizeType ind, String name)
	{
		if (ind >= value_.size())
			throw RuntimeError("No value for " + name + "\n");

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

	void setVectorValue(VectorStringType& v, String rhs, String name)
	{
		SizeType last = rhs.length();
		if (last < 2)
			err("Vector must be enclosed in brakets, name= " + name + "\n");

		--last;
		if (rhs[0] != '[' || rhs[last] != ']')
			err("Vector must be enclosed in brakets, name= " + name + "\n");

		rhs = (last < 2 ) ? "" : rhs.substr(1,last - 1);
		AinurLexicalType::removeTrailingBlanks(rhs);
		last = rhs.length();
		if (last > 1 && rhs[0] == '[' && rhs[--last] == ']') {
			// it's really a matrix
			setMatrixValue("[" + rhs + "]", name);
			type_ = MATRIX;
			return;
		}

		split(v, rhs, ",");
	}

	// [[a, b, c], [a, b, c]]
	void setMatrixValue(String rhs, String name)
	{
		SizeType last = rhs.length();
		if (last < 4)
			err("Matrix must be enclosed in brakets\n");

		--last;
		if (rhs[0] != '[' || rhs[last] != ']')
			err("Matrix must be enclosed in brakets " + rhs + "\n");
		assert(last > 2);
		rhs =  rhs.substr(1,last - 1);
		VectorStringType tmp;
		split(tmp, rhs, "[");
		// a, b, c],
		// a, b, c]
		SizeType rows = tmp.size();
		assert(rows > 0);
		SizeType cols = 0;
		value_.clear();
		SizeType offset = 2;
		for (SizeType row = 0; row < rows; ++row) {
			VectorStringType v;
			String s = tmp[row];
			AinurLexicalType::removeTrailingBlanks(s);
			SizeType last = s.length();
			if (last > 0 && s[--last] == ',')
				s = s.substr(0, last);

			s = "[" + s;
			setVectorValue(v, s, name);
			SizeType thisCol = v.size();
			if (row == 0) {
				cols = thisCol;
				value_.resize(rows*cols + 2);
				value_[0] = ttos(rows);
				value_[1] = ttos(cols);
			} else if (cols != thisCol) {
				err("Malformed matrix, " + name + "\n");
			}

			appendToVecStr(value_, v, offset);
			offset += v.size();
		}
	}

	void appendToVecStr(VectorStringType& dest,
	                    const VectorStringType& src,
	                    SizeType offset) const
	{
		SizeType n = src.size();
		assert(offset + n <= dest.size());
		for (SizeType i = 0; i < n; ++i)
			dest[offset + i] = src[i];
	}

	Type type_;
	SubType subType_;
	Attribute attr_;
	VectorStringType value_;
	mutable SizeType used_;
};
}
#endif // AINURSTORE_H
