#ifndef STRINGORDERPOST_H
#define STRINGORDERPOST_H
#include "Vector.h"
#include "Braket.h"

namespace Dmrg {

template<typename ModelType>
class StringOrderPost {

public:

	typedef Braket<ModelType> BraketType;
	typedef typename BraketType::VectorStringType VectorStringType;
	typedef typename BraketType::AlgebraType AlgebraType;
	typedef typename BraketType::OperatorSpecType OperatorSpecType;
	typedef typename BraketType::VectorAlgebraType VectorAlgebraType;
	typedef typename BraketType::VectorIntType VectorIntType;
	typedef typename BraketType::MatrixType MatrixType;

	StringOrderPost(const BraketType& braket)
	{
		static const PsimagLite::String stringop = "!stringorder=";
		const PsimagLite::String special = braket.opName(0);

		const SizeType l = stringop.length();
		if (special.substr(0, l) == stringop) {
			PsimagLite::String rest = special.substr(l, special.length() - l);
			VectorStringType vStr;
			PsimagLite::split(vStr, rest, ":");
			addAllOps(vStr, stringop, braket.model());
			return;
		}

		err("FATAL: Braket: special op " + special + " not understood\n");
	}

	void computeMatrix(MatrixType& m)
	{
		err("testing\n");
	}

private:

	void addAllOps(const VectorStringType& vStr,
	               PsimagLite::String stringop,
	               const ModelType& model)
	{
		if (vStr.size() != 3)
			err("FATAL: Braket: " + stringop + " expecting 3 elements to follow\n");

		VectorIntType sites(vStr.size(), -1);
		OperatorSpecType opSpec(model);
		PsimagLite::CanonicalExpression<OperatorSpecType> canonicalExpression(opSpec);
		const AlgebraType opEmpty;
		for (SizeType i = 0; i < vStr.size(); ++i) {
			AlgebraType tmp;
			canonicalExpression(tmp, vStr[i], opEmpty, sites[i]);
			op_.push_back(tmp);
		}
	}

	VectorAlgebraType op_;
};
}
#endif // STRINGORDERPOST_H
