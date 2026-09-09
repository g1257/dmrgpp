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

namespace PsimagLite {

template <typename VectorValueType> class Modulus : public Node<VectorValueType> {

	using ValueType = typename VectorValueType::value_type;

public:

	Modulus* clone() const override { return new Modulus(*this); }

	PsimagLite::String code() const override { return "%"; }

	SizeType arity() const override { return 2; }

	ValueType exec(const VectorValueType& v) const override
	{
		using RealType = typename Real<ValueType>::Type;

		assert(v.size() == 2);

		RealType v0 = PsimagLite::norm(v[0]);
		RealType v1 = PsimagLite::norm(v[1]);

		SizeType x = static_cast<int>(v0);
		SizeType y = static_cast<int>(v1);

		ValueType result = (y == 0) ? x : x % y;
		return result;
	}

}; // class Modulus

template <typename VectorValueType> class Cosine : public Node<VectorValueType> {

	using ValueType = typename VectorValueType::value_type;

public:

	Cosine* clone() const override { return new Cosine(*this); }

	PsimagLite::String code() const override { return "cos"; }

	SizeType arity() const override { return 1; }

	ValueType exec(const VectorValueType& v) const override
	{
		assert(v.size() == 1);
		return cos(v[0]);
	}
}; // class Cosine

template <typename VectorValueType> class Sine : public Node<VectorValueType> {

	using ValueType = typename VectorValueType::value_type;

public:

	Sine* clone() const override { return new Sine(*this); }

	PsimagLite::String code() const override { return "sin"; }

	SizeType arity() const override { return 1; }

	ValueType exec(const VectorValueType& v) const override
	{
		assert(v.size() == 1);
		return sin(v[0]);
	}
}; // class Sine

template <typename VectorValueType> class Exp : public Node<VectorValueType> {

	using ValueType = typename VectorValueType::value_type;

public:

	Exp* clone() const override { return new Exp(*this); }

	PsimagLite::String code() const override { return "exp"; }

	SizeType arity() const override { return 1; }

	ValueType exec(const VectorValueType& v) const override
	{
		assert(v.size() == 1);
		return std::exp(v[0]);
	}
}; // class Exp

template <typename VectorValueType> class TernaryOp : public Node<VectorValueType> {

	using ValueType = typename VectorValueType::value_type;

public:

	TernaryOp* clone() const override { return new TernaryOp(*this); }

	PsimagLite::String code() const override { return "?"; }

	SizeType arity() const override { return 3; }

	ValueType exec(const VectorValueType& v) const override
	{
		assert(v.size() == 3);
		SizeType b = static_cast<SizeType>(PsimagLite::norm(v[0]));
		return (b) ? v[1] : v[2];
	}
}; // class TernaryOp

template <typename VectorValueType> class Log : public Node<VectorValueType> {

	using ValueType = typename VectorValueType::value_type;

public:

	Log* clone() const override { return new Log(*this); }

	PsimagLite::String code() const override { return "log"; }

	SizeType arity() const override { return 1; }

	ValueType exec(const VectorValueType& v) const override
	{
		assert(v.size() == 1);
		return log(v[0]);
	}
}; // class Log
} // namespace PsimagLite

#endif // ADDITIONALFUNCTIONS_H
