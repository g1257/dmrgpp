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
#ifndef PLUS_MINUS_MULT_DIV_H
#define PLUS_MINUS_MULT_DIV_H
#include "../Vector.h"
#include <cassert>
#include "Node.h"
#include "../PsimagLite.h"
#include "AdditionalFunctions.h"

namespace PsimagLite {

template<typename ValueType_>
class PlusMinusMultiplyDivide {

public:

	typedef typename PsimagLite::Vector<ValueType_>::Type VectorValueType;
	typedef Node<VectorValueType> NodeType;
	typedef typename PsimagLite::Vector<NodeType*>::Type VectorNodeType;
	typedef Plus<VectorValueType> PlusType;
	typedef Minus<VectorValueType> MinusType;
	typedef Times<VectorValueType> TimesType;
	typedef DividedBy<VectorValueType> DividedByType;
	typedef Modulus<VectorValueType> ModulusType;
	typedef Input<VectorValueType> InputType;
	typedef ValueType_ ValueType;
	typedef PsimagLite::Vector<PsimagLite::String>::Type VectorStringType;
	typedef PsimagLite::Vector<SizeType>::Type VectorSizeType;

	PlusMinusMultiplyDivide() : inputsSoFar_(0)
	{
		NodeType* plus = new PlusType();
		nodes_.push_back(plus);

		NodeType* minus = new MinusType();
		nodes_.push_back(minus);

		NodeType* times = new TimesType();
		nodes_.push_back(times);

		NodeType* dividedBy = new DividedByType();
		nodes_.push_back(dividedBy);

		NodeType* modulus = new ModulusType();
		nodes_.push_back(modulus);

		NodeType* sine = new Sine<VectorValueType>();
		nodes_.push_back(sine);

		NodeType* cosine = new Cosine<VectorValueType>();
		nodes_.push_back(cosine);

		NodeType* log = new Log<VectorValueType>();
		nodes_.push_back(log);

		NodeType* exp = new Exp<VectorValueType>();
		nodes_.push_back(exp);

		NodeType* ternary = new TernaryOp<VectorValueType>();
		nodes_.push_back(ternary);
	}

	~PlusMinusMultiplyDivide()
	{
		for (SizeType i = 0; i < nodes_.size(); i++)
			delete nodes_[i];

		nodes_.clear();
	}

	const VectorNodeType& nodesSerial() const
	{
		return nodes_;
	}

	const NodeType& findNodeFromCode(const String& code)
	{
		SizeType n = nodes_.size();
		for (SizeType i = 0; i < n; ++i) {
			if (nodes_[i]->isInput()) continue;
			if (nodes_[i]->code() == code) return *nodes_[i];
		}

		// Assume it's an input
		NodeType* input = new InputType(inputsSoFar_++);

		try {
		input->set(PsimagLite::atof(code));
		} catch (std::exception&) {
			throw RuntimeError("FATAL: No function with code " + code + " was found.\n");
		}

		nodes_.push_back(input);
		return *input;
	}

private:

	VectorNodeType nodes_;
	SizeType inputsSoFar_;
}; // class PlusMinusMultiplyDivide

} // namespace PsimagLite

#endif // PLUS_MINUS_MULT_DIV_H
