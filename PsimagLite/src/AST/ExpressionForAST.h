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
#ifndef EXPRESSIONFORAST_H
#define EXPRESSIONFORAST_H
#include "../Vector.h"
#include "Node.h"
#include "Tree.h"

namespace PsimagLite {

template <typename PrimitivesType>
class ExpressionForAST {
public:

	typedef typename PrimitivesType::VectorValueType VectorValueType;
	typedef typename VectorValueType::value_type ValueType;
	typedef Tree<PrimitivesType> TreeType;
	typedef Vector<String>::Type VectorStringType;
	typedef Node<VectorValueType> NodeType;
	typedef typename Vector<TreeType*>::Type VectorTreeType;

	ExpressionForAST(const VectorStringType& vecStr,
	                 PrimitivesType& primitives)
	{
		constexpr bool verbose = false; // FIXME
		SizeType effectiveSize = vecStr.size();
		PsimagLite::Vector<SizeType>::Type va;

		SizeType sumOfA = 1;
		for (SizeType i = 0; i < effectiveSize; i++) {
			PsimagLite::String cStr = vecStr[i];
			const NodeType& node = primitives.findNodeFromCode(cStr);

			SizeType a = node.arity();
			sumOfA += (a - 1);
			TreeType* tree = new TreeType(primitives, node, verbose);

			va.push_back(a);
			vecTree_.push_back(tree);
			if (sumOfA == 0)
				break;
		}

		SizeType k = 0;
		for (SizeType i = 0; i < vecTree_.size(); i++) {
			SizeType a = va[i];
			if (a == 0 || !vecTree_[i])
				continue;
			for (SizeType j = k + 1; j < k + a + 1; j++) {
				if (j >= vecTree_.size())
					continue;
				vecTree_[i]->setDescendants(*vecTree_[j]);
			}

			k += a;
		}
	}

	~ExpressionForAST()
	{
		for (SizeType i = 0; i < vecTree_.size(); i++)
			delete vecTree_[i];

		vecTree_.clear();
	}

	ValueType exec()
	{
		assert(0 < vecTree_.size());
		assert(vecTree_[0]);
		return vecTree_[0]->exec();
	}

private:

	VectorTreeType vecTree_;
}; // class ExpressionForAST

} // namespace PsimagLite
#endif // EXPRESSIONFORAST_H
