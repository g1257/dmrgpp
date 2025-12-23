/*
Copyright (c) 2009-2022, UT-Battelle, LLC
All rights reserved

[PsimagLite, Version 2.]

*********************************************************
THE SOFTWARE IS SUPPLIED BY THE COPYRIGHT HOLDERS AND
CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED
WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A
PARTICULAR PURPOSE ARE DISCLAIMED.

Please see full open source license included in file LICENSE.
*********************************************************
*/
#ifndef PSI_NODE_H
#define PSI_NODE_H
#include "../Vector.h"
#include <cassert>

namespace PsimagLite {

template <typename VectorValueType, typename AnglesType_ = int> class Node {

public:

	typedef AnglesType_ AnglesType;
	typedef typename VectorValueType::value_type ValueType;
	typedef typename PsimagLite::Vector<AnglesType>::Type VectorAnglesType;

	virtual ~Node() { }

	virtual Node* clone() const = 0;

	virtual PsimagLite::String code() const = 0;

	virtual SizeType arity() const = 0;

	virtual ValueType exec(const VectorValueType&) const = 0;

	virtual ValueType exec(const VectorValueType&, const VectorAnglesType*, SizeType&) const
	{
		throw PsimagLite::RuntimeError("node::exec() long form\n");
	}

	virtual void set(const ValueType&) const { throw PsimagLite::RuntimeError("node::set\n"); }

	virtual void setAngle(PsimagLite::String) const { }

	virtual void print(std::ostream&) const { }

	virtual void setDcValue(const ValueType&) const
	{
		throw PsimagLite::RuntimeError("node::setDcValue\n");
	}

	virtual bool isInput() const { return false; }

}; // class Node

template <typename VectorValueType> class Plus : public Node<VectorValueType> {

	typedef typename VectorValueType::value_type ValueType;

public:

	Plus* clone() const { return new Plus(*this); }

	virtual PsimagLite::String code() const { return "+"; }

	virtual SizeType arity() const { return 2; }

	virtual ValueType exec(const VectorValueType& v) const
	{
		assert(v.size() == 2);
		return v[0] + v[1];
	}

}; // class Plus

template <typename VectorValueType> class Minus : public Node<VectorValueType> {

	typedef typename VectorValueType::value_type ValueType;

public:

	Minus* clone() const { return new Minus(*this); }

	virtual PsimagLite::String code() const { return "-"; }

	virtual SizeType arity() const { return 2; }

	virtual ValueType exec(const VectorValueType& v) const
	{
		assert(v.size() == 2);
		return v[0] - v[1];
	}

}; // class Minus

template <typename VectorValueType> class Times : public Node<VectorValueType> {

	typedef typename VectorValueType::value_type ValueType;

public:

	Times* clone() const { return new Times(*this); }

	virtual PsimagLite::String code() const { return "*"; }

	virtual SizeType arity() const { return 2; }

	virtual ValueType exec(const VectorValueType& v) const
	{
		assert(v.size() == 2);
		return v[0] * v[1];
	}

}; // class Times

template <typename VectorValueType> class DividedBy : public Node<VectorValueType> {

	typedef typename VectorValueType::value_type ValueType;

public:

	DividedBy* clone() const { return new DividedBy(*this); }

	virtual PsimagLite::String code() const { return "/"; }

	virtual SizeType arity() const { return 2; }

	virtual ValueType exec(const VectorValueType& v) const
	{
		assert(v.size() == 2);
		if (std::norm(v[1]) < 1e-6)
			return v[0];

		return v[0] / v[1];
	}

}; // class DividedBy

template <typename VectorValueType> class IfGtZero : public Node<VectorValueType> {

	typedef typename VectorValueType::value_type ValueType;

public:

	IfGtZero* clone() const { return new IfGtZero(*this); }

	virtual char code() const { return 'g'; }

	virtual SizeType arity() const { return 1; }

	virtual ValueType exec(const VectorValueType& v) const
	{
		assert(v.size() == 1);

		return (v[0] > 0) ? 1 : 0;
	}

}; // class IfGtZero

template <typename VectorValueType> class Int : public Node<VectorValueType> {

	typedef typename VectorValueType::value_type ValueType;

public:

	Int* clone() const { return new Int(*this); }

	virtual char code() const { return 'i'; }

	virtual SizeType arity() const { return 1; }

	virtual ValueType exec(const VectorValueType& v) const
	{
		assert(v.size() == 1);

		SizeType x = static_cast<SizeType>(v[0]);
		return x;
	}

}; // class Int

template <typename VectorValueType> class Input : public Node<VectorValueType> {

	typedef typename VectorValueType::value_type ValueType;

public:

	Input(SizeType i)
	    : char_(i + 48)
	    , strOneChar_(" ")
	{
		strOneChar_[0] = char_;
	}

	Input* clone() const { return new Input(*this); }

	virtual PsimagLite::String code() const { return strOneChar_; }

	virtual SizeType arity() const { return 0; }

	virtual ValueType exec(const VectorValueType& v) const
	{
		assert(v.size() == 0);
		return input_;
	}

	virtual void set(const ValueType& x) const { input_ = x; }

	virtual bool isInput() const { return true; }

private:

	char char_;
	PsimagLite::String strOneChar_;
	mutable ValueType input_;

}; // class Input

} // namespace PsimagLite
#endif // PSI_NODE_H
