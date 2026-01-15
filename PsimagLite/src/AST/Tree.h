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

#ifndef PSI_TREE_H
#define PSI_TREE_H
#include "../PsimagLite.h"

namespace PsimagLite {

template <typename PrimitivesType> class Tree {

public:

	typedef Tree<PrimitivesType>                               TreeType;
	typedef typename PrimitivesType::ValueType                 ValueType;
	typedef typename PrimitivesType::NodeType                  NodeType;
	typedef typename PsimagLite::Vector<const TreeType*>::Type VectorTreeType;
	typedef typename NodeType::AnglesType                      AnglesType;
	typedef typename PsimagLite::Vector<AnglesType>::Type      VectorAnglesType;
	typedef typename PsimagLite::Vector<ValueType>::Type       VectorValueType;

	Tree(const PrimitivesType& primitives, const NodeType& node, bool verbose)
	    : primitives_(primitives)
	    , node_(node)
	    , verbose_(verbose)
	{ }

	~Tree() { }

	ValueType exec() const
	{
		if (verbose_)
			std::cout << " type= " << node_.code() << "\n";
		VectorValueType values(descendants_.size());

		for (SizeType i = 0; i < descendants_.size(); i++) {
			values[i] = descendants_[i]->exec();
		}

		//		for (SizeType i = 0; i < values.size(); i++) {
		//			if (values[i]<0 || values[i]==0 ||
		// values[i] >0) continue; 			throw
		// std::runtime_error("exec\n");
		//		}

		ValueType tmp = node_.exec(values);
		if (verbose_) {
			std::cout << "tmp= " << tmp << " type= " << node_.code();
			for (SizeType i = 0; i < values.size(); i++)
				std::cout << values[i] << " ";
			std::cout << "\n";
		}

		return tmp;
	}

	void set(const VectorValueType& values) const
	{
		for (SizeType i = 0; i < descendants_.size(); i++)
			descendants_[i]->set(values);

		const PsimagLite::String str = node_.code();
		if (str.length() == 0)
			err("Node with empty code!?\n");

		if (str.length() > 1)
			return;
		const char c = str[0];
		if (c < 48 || c > 57)
			return;
		SizeType index = c - 48;
		assert(index < values.size());
		node_.set(values[index]);
	}

	bool isLinearTree() const
	{
		if (descendants_.size() > 1)
			return false;

		if (descendants_.size() == 0)
			return true;

		return descendants_[0]->isLinearTree();
	}

	void setDescendants(const TreeType& n0) { descendants_.push_back(&n0); }

	void setDescendants(const TreeType& n0, const TreeType& n1)
	{
		descendants_.push_back(&n0);
		descendants_.push_back(&n1);
	}

private:

	Tree& operator=(const Tree&) = delete;

	const PrimitivesType& primitives_;
	const NodeType&       node_;
	bool                  verbose_;
	VectorTreeType        descendants_;
};

} // namespace PsimagLite

#endif // PSI_TREE_H
