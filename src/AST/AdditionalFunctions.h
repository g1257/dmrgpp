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

#ifndef ADDITIONALFUNCTIONS_H
#define ADDITIONALFUNCTIONS_H
#include "Node.h"

namespace PsimagLite
{

template <typename VectorValueType>
class Modulus : public Node<VectorValueType>
{

	typedef typename VectorValueType::value_type ValueType;

public:

	Modulus* clone() const { return new Modulus(*this); }

	virtual PsimagLite::String code() const { return "%"; }

	virtual SizeType arity() const { return 2; }

	virtual ValueType exec(const VectorValueType& v) const
	{
		typedef typename Real<ValueType>::Type RealType;

		assert(v.size() == 2);

		RealType v0 = PsimagLite::norm(v[0]);
		RealType v1 = PsimagLite::norm(v[1]);

		SizeType x = static_cast<int>(v0);
		SizeType y = static_cast<int>(v1);

		ValueType result = (y == 0) ? x : x % y;
		return result;
	}

}; // class Modulus

template <typename VectorValueType>
class Cosine : public Node<VectorValueType>
{

	typedef typename VectorValueType::value_type ValueType;

public:

	Cosine* clone() const { return new Cosine(*this); }

	virtual PsimagLite::String code() const { return "cos"; }

	virtual SizeType arity() const { return 1; }

	virtual ValueType exec(const VectorValueType& v) const
	{
		assert(v.size() == 1);
		return cos(v[0]);
	}
}; // class Cosine

template <typename VectorValueType>
class Sine : public Node<VectorValueType>
{

	typedef typename VectorValueType::value_type ValueType;

public:

	Sine* clone() const { return new Sine(*this); }

	virtual PsimagLite::String code() const { return "sin"; }

	virtual SizeType arity() const { return 1; }

	virtual ValueType exec(const VectorValueType& v) const
	{
		assert(v.size() == 1);
		return sin(v[0]);
	}
}; // class Sine

template <typename VectorValueType>
class Exp : public Node<VectorValueType>
{

	typedef typename VectorValueType::value_type ValueType;

public:

	Exp* clone() const { return new Exp(*this); }

	virtual PsimagLite::String code() const { return "exp"; }

	virtual SizeType arity() const { return 1; }

	virtual ValueType exec(const VectorValueType& v) const
	{
		assert(v.size() == 1);
		return std::exp(v[0]);
	}
}; // class Exp

template <typename VectorValueType>
class TernaryOp : public Node<VectorValueType>
{

	typedef typename VectorValueType::value_type ValueType;

public:

	TernaryOp* clone() const { return new TernaryOp(*this); }

	virtual PsimagLite::String code() const { return "?"; }

	virtual SizeType arity() const { return 3; }

	virtual ValueType exec(const VectorValueType& v) const
	{
		assert(v.size() == 3);
		SizeType b = static_cast<SizeType>(PsimagLite::norm(v[0]));
		return (b) ? v[1] : v[2];
	}
}; // class TernaryOp

template <typename VectorValueType>
class Log : public Node<VectorValueType>
{

	typedef typename VectorValueType::value_type ValueType;

public:

	Log* clone() const { return new Log(*this); }

	virtual PsimagLite::String code() const { return "log"; }

	virtual SizeType arity() const { return 1; }

	virtual ValueType exec(const VectorValueType& v) const
	{
		assert(v.size() == 1);
		return log(v[0]);
	}
}; // class Log
} // namespace PsimagLite

#endif // ADDITIONALFUNCTIONS_H
