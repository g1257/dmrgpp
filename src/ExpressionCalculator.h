#ifndef PSI_EXPRESSIONCALCULATOR_H
#define PSI_EXPRESSIONCALCULATOR_H

#include <iostream>
#include <stack>
#include "Vector.h"
#include "Tokenizer.h"
#include "TypeToString.h"
#include "Stack.h"

namespace PsimagLite {

template<typename ComplexOrRealType>
struct PrepassData {
	typedef typename PsimagLite::Vector<ComplexOrRealType>::Type VectorType;

	PsimagLite::String names;
	VectorType values;
};

template<typename PrepassDataType>
class ExpressionPrepass {

	typedef PsimagLite::Vector<PsimagLite::String>::Type VectorStringType;

public:

	static void prepass(VectorStringType& vs,
	                    const PrepassDataType& pd)
	{
		for (SizeType i = 0; i < vs.size(); ++i) {
			int ri = getReplacementIndex(vs[i],pd);
			if (ri < 0) continue;
			vs[i] = ttos(pd.values[ri]);
		}
	}

private:

	static int getReplacementIndex(PsimagLite::String str,
	                               const PrepassDataType& pd)
	{
		SizeType l = str.size();
		if (l < 2) return -1;
		if (str[0] != '%') return -1;
		unsigned char letter = str[1];
		for (SizeType i = 0; i < pd.names.size(); ++i)
			if (pd.names[i] == letter) return i;
		return -1;
	}
}; // class ExpressionPrepass

template<typename ComplexOrRealType>
class ExpressionCalculator {

	typedef PsimagLite::Vector<SizeType>::Type VectorSizeType;
	typedef typename PsimagLite::Vector<ComplexOrRealType>::Type VectorType;

	enum NodeTypeEnum {NODE_NUMBER, NODE_OPERATOR};

	static const bool debug_ = false;

	struct Node {

		Node(PsimagLite::String str,SizeType ind)
		: type(NODE_OPERATOR),index(ind),ary(0),value(0.0)
		{
			for (SizeType i = 0; i < 4; ++i) op[i] = '\0';
			SizeType l = str.size();
			if (l == 0) return;
			if (str[l-1] >= '0' && str[l-1]<='9') {
				type = NODE_NUMBER;
				value = atof(str.c_str());
			} else {
				ary = findAry(str);
				for (SizeType i = 0; i < 4; ++i)
					op[i] = (i < l) ? str[i] : '\0';
			}
		}

		NodeTypeEnum type;
		unsigned char op[4];
		SizeType index;
		SizeType ary;
		ComplexOrRealType value;

		void print(std::ostream& os) const
		{
			if (type == NODE_NUMBER) {
				os<<index<<" "<<value<<"| ";
			} else {
				os<<index<<" "<<op<<" "<<ary<<"| ";
			}
		}
	};

	typedef Node NodeType;
	typedef typename PsimagLite::Vector<NodeType>::Type VectorNodeType;

public:

	typedef PsimagLite::Vector<PsimagLite::String>::Type VectorStringType;

	ExpressionCalculator(const VectorStringType& ve)
	: value_(0)
	{
		VectorNodeType vne(ve.size(),NodeType("",0));
		fillNodes(vne,ve);
		while (!onePass(vne)) {};
	}

	bool onePass(VectorNodeType& vne)
	{
		typename Stack<Node>::Type ops;
		SizeType numbers = 0;
		VectorSizeType indices(3,0);
		SizeType simplifications = 0;
		for (SizeType i = 0; i < vne.size(); ++i) {
			if (vne[i].type == NODE_NUMBER) {
				if (ops.size() == 0) {
					value_ = vne[i].value;
					return true;
				} else {
					indices[numbers] = i;
					numbers++;
					if (ops.top().ary == numbers) {
						simplify(vne,ops.top(),indices);
						numbers = 0;
						ops.pop();
						simplifications++;
					}
				}
			} else {
				if (vne[i].op[0] != '\0')
					ops.push(vne[i]);
			}
		}

		if (debug_) print(vne,"vne");
		if (simplifications == 0)
			throw PsimagLite::RuntimeError("Syntax Error\n");

		return false;
	}

	const ComplexOrRealType& operator()() const
	{
		return value_;
	}

private:

	void print(const VectorNodeType& vn, PsimagLite::String label) const
	{
		std::cout<<label<<" ";
		for (SizeType i = 0; i < vn.size(); ++i)
			vn[i].print(std::cout);
		std::cout<<"\n";
	}

	void fillNodes(VectorNodeType& vn, const VectorStringType& ve) const
	{
		for (SizeType i = 0; i < ve.size(); ++i)
			vn[i] = Node(ve[i],i);
	}

	void simplify(VectorNodeType& vn,
	              const NodeType& op,
	              const VectorSizeType& indices) const
	{
		SizeType total = op.ary;
		VectorType values(total);
		for (SizeType i = 0; i < total; ++i) {
			values[i] = vn[indices[i]].value;
			vn[indices[i]] = Node("",indices[i]);
		}

		ComplexOrRealType value = executeOperator(op.op,values);
		vn[op.index] = Node(ttos(value),op.index);
	}

	static ComplexOrRealType executeOperator(const unsigned char op[],
	                                const VectorType& values)
	{
		if (op[0] == '+') return values[0] + values[1];
		if (op[0] == '-') return values[0] - values[1];
		if (op[0] == '*') return values[0] * values[1];
		if (op[0] == 'c') return cos(values[0]);
		if (op[0] == 's') return sin(values[0]);
		if (op[0] == '?') return (std::real(values[0]) > 0) ? values[1] : values[2];
		return 0.0;
	}

	static SizeType findAry(PsimagLite::String op)
	{
		if (op == "+" || op == "-" || op == "*") return 2;
		if (op == "c" || op == "s") return 1;
		if (op == "?") return 3;
		return 0;
	}

	ComplexOrRealType value_;
}; // class ExpressionCalculator

} // namespace PsimagLite
#endif // PSI_EXPRESSIONCALCULATOR_H

