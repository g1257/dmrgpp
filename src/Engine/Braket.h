#ifndef DMRG_braket_H
#define DMRG_braket_H
#include "CanonicalExpression.h"
#include "GetBraOrKet.h"
#include "Matrix.h"
#include "OperatorSpec.h"
#include "PsimagLite.h"
#include "Vector.h"

namespace Dmrg
{

template <typename ModelType_>
class Braket
{

public:

	typedef ModelType_ ModelType;
	typedef typename ModelType::OperatorType OperatorType;
	typedef typename PsimagLite::Vector<int>::Type VectorIntType;
	typedef typename OperatorType::PairType PairType;
	typedef typename OperatorType::Su2RelatedType Su2RelatedType;
	typedef typename OperatorType::StorageType SparseMatrixType;
	typedef typename SparseMatrixType::value_type ComplexOrRealType;
	typedef PsimagLite::Matrix<ComplexOrRealType> MatrixType;
	typedef typename PsimagLite::Real<ComplexOrRealType>::Type RealType;
	typedef OperatorSpec<ModelType, OperatorType> OperatorSpecType;
	typedef typename OperatorSpecType::ResultType AlgebraType;
	typedef typename OperatorSpecType::OneOperatorSpecType OneOperatorSpecType;
	typedef typename PsimagLite::Vector<AlgebraType>::Type VectorAlgebraType;
	typedef PsimagLite::Vector<PsimagLite::String>::Type VectorStringType;
	typedef PsimagLite::GetBraOrKet GetBraOrKetType;
	typedef PsimagLite::Vector<GetBraOrKetType>::Type VectorGetBraOrKetType;

	Braket(const ModelType& model, const PsimagLite::String& braket)
	    : model_(model)
	    , savedString_(braket)
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

		braket_.push_back(vecStr[0].substr(1, vecStr[0].length() - 1));

		braket_.push_back(vecStr[2].substr(0, vecStr[2].length() - 1));

		if (vecStr[1].length() > 1 && vecStr[1][0] == '!') {
			opExprName_.resize(1);
			opExprName_[0] = vecStr[1];
			sites_.resize(opExprName_.size(), -1);
			return; // early exit <===
		}

		PsimagLite::split(opExprName_, vecStr[1], ";");

		sites_.resize(opExprName_.size(), -1);

		OperatorSpecType opSpec(model);
		PsimagLite::CanonicalExpression<OperatorSpecType> canonicalExpression(opSpec);
		const AlgebraType opEmpty;
		for (SizeType i = 0; i < opExprName_.size(); ++i) {
			AlgebraType tmp;
			canonicalExpression(tmp, opExprName_[i], opEmpty, sites_[i]);
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

	const GetBraOrKetType& bra() const
	{
		assert(braket_.size() > 0);
		return braket_[0];
	}

	const GetBraOrKetType& ket() const
	{
		assert(braket_.size() > 1);
		return braket_[1];
	}

	SizeType points() const { return opExprName_.size(); }

	SizeType site(SizeType ind) const
	{
		assert(ind < sites_.size());
		if (sites_[ind] >= 0)
			return sites_[ind];
		throw PsimagLite::RuntimeError("site is negative\n");
	}

	PsimagLite::String toString() const { return savedString_; }

	const ModelType& model() const { return model_; }

	// avoid using this function in new code
	// it's only to support legacy code
	void forceOperators(const OperatorType& op1, const OperatorType& op2)
	{
		op_.resize(2);
		op_[0] = op1;
		op_[1] = op2;
	}

	void reorder(const std::vector<SizeType>& permutation)
	{
		SizeType n = permutation.size();
		if (n != op_.size()) {
			err("Braket::reorder() size mismatch\n");
		}

		if (op_.size() != sites_.size()) {
			err("Braket::reorder() size mismatch; consistency error\n");
		}

		if (n == 1)
			return;

		std::string str = "<" + braket_[0].toString() + "|";
		VectorAlgebraType opNew(n);
		VectorIntType sitesNew(n);
		VectorStringType opExprNameNew(n);
		for (SizeType i = 0; i < n; ++i) {
			opNew[i] = op_[permutation[i]];
			sitesNew[i] = sites_[permutation[i]];
			opExprNameNew[i] = opExprName_[permutation[i]];
			str += opName(i);
			if (i + 1 < n)
				str += std::string(";");
		}

		str += braket_[1].toString() + ">";

		// actual changes
		savedString_ = str;
		op_ = opNew;
		sites_ = sitesNew;
		opExprName_ = opExprNameNew;
		// braket_ does not change because bra and ket don't change
	}

private:

	const ModelType& model_;
	VectorGetBraOrKetType braket_; // convert to pair?
	PsimagLite::String savedString_;
	VectorStringType opExprName_;
	VectorAlgebraType op_;
	VectorIntType sites_;
}; // class Braket
}

#endif // DMRG_braket_H
