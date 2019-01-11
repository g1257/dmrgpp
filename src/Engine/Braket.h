#ifndef DMRG_braket_H
#define DMRG_braket_H
#include "Vector.h"
#include "PsimagLite.h"
#include "Matrix.h"
#include "OperatorSpec.h"
#include "CanonicalExpression.h"

namespace Dmrg {

template<typename ModelType>
class Braket {

	typedef typename ModelType::OperatorType OperatorType;
	typedef typename PsimagLite::Vector<int>::Type VectorIntType;
	typedef typename OperatorType::PairType PairType;
	typedef typename OperatorType::Su2RelatedType Su2RelatedType;
	typedef typename OperatorType::StorageType SparseMatrixType;
	typedef typename SparseMatrixType::value_type ComplexOrRealType;
	typedef PsimagLite::Matrix<ComplexOrRealType> MatrixType;
	typedef typename PsimagLite::Real<ComplexOrRealType>::Type RealType;
	typedef OperatorSpec<ModelType, OperatorType> OperatorSpecType;

public:

	typedef typename OperatorSpecType::ResultType AlgebraType;
	typedef typename PsimagLite::Vector<AlgebraType>::Type VectorAlgebraType;
	typedef PsimagLite::Vector<PsimagLite::String>::Type VectorStringType;

	Braket(const ModelType& model,const PsimagLite::String& braket)
	    : model_(model), braket_(2,""),savedString_(braket)
	{
		VectorStringType vecStr;
		PsimagLite::split(vecStr, braket, "|");

		if (vecStr.size() == 2) {
			VectorStringType tmp;
			tmp.push_back(vecStr[0]);
			tmp.push_back("identity");
			tmp.push_back(vecStr[1]);
			vecStr = tmp;
		}

		if (vecStr.size() != 3) {
			PsimagLite::String str("ObserverInterpreter: syntax error for ");
			str += braket + " is not a Braket\n";
			throw PsimagLite::RuntimeError(str);
		}

		braket_[0] = vecStr[0].substr(1,vecStr[0].length()-1);
		if (!isBraket(0)) {
			PsimagLite::String str("ObserverInterpreter: syntax error: ");
			str += braket_[0] + " must be <gs or <time or <P\\d+\n";
			throw PsimagLite::RuntimeError(str);
		}

		braket_[1] = vecStr[2].substr(0,vecStr[2].length()-1);
		if (!isBraket(1)) {
			PsimagLite::String str("ObserverInterpreter: syntax error: ");
			str += braket_[1] + " must be gs> or time>  or <P\\d+\n";
			throw PsimagLite::RuntimeError(str);
		}

		PsimagLite::split(opExprName_, vecStr[1], ";");

		sites_.resize(opExprName_.size(),-1);

		OperatorSpecType opSpec(model);
		PsimagLite::CanonicalExpression<OperatorSpecType> canonicalExpression(opSpec);

		for (SizeType i = 0; i < opExprName_.size(); ++i) {
			AlgebraType tmp = canonicalExpression(opExprName_[i],sites_[i]);
			op_.push_back(tmp);
		}
	}

	const AlgebraType& op(SizeType ind) const
	{
		assert(ind < op_.size());
		return op_[ind];
	}

	PsimagLite::String opName(SizeType ind) const
	{
		assert(ind < opExprName_.size());
		return opExprName_[ind];
	}

	PsimagLite::String bra() const
	{
		assert(isBraket(0));
		return braket_[0];
	}

	PsimagLite::String ket() const
	{
		assert(isBraket(1));
		return braket_[1];
	}

	SizeType points() const { return opExprName_.size(); }

	SizeType site(SizeType ind) const
	{
		if (sites_[ind] >= 0) return sites_[ind];
		throw PsimagLite::RuntimeError("site is negative\n");
	}

	PsimagLite::String toString() const { return savedString_; }

	static int getPtype(PsimagLite::String str)
	{
		// str == P\d+
		if (str.length() < 2) return -1;
		if (str[0] != 'P') return -1;
		PsimagLite::String number("");
		for (SizeType i = 1; i < str.length(); ++i) {
			number += str[i];
			unsigned char x = str[i];
			if (x < 48 || x > 57) return -1;
		}

		return atoi(number.c_str()) + 1;
	}

private:

	bool isBraket(SizeType ind) const
	{
		if (ind >= braket_.size()) return false;
		int pType = getPtype(braket_[ind]);

		return (braket_[ind] == "gs" ||
		        braket_[ind] == "time" ||
		        pType >= 0);
	}

	const ModelType& model_;
	VectorStringType braket_;
	PsimagLite::String savedString_;
	VectorStringType opExprName_;
	SizeType type_;
	VectorAlgebraType op_;
	VectorIntType sites_;
}; // class Braket
}

#endif // DMRG_braket_H

