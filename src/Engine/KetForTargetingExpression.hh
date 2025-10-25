#ifndef KETFORTARGETINGEXPRESSION_HH
#define KETFORTARGETINGEXPRESSION_HH
#include "PsimagLite.h"
#include "FactorForTargetingExpression.hh"

namespace Dmrg {

template<typename ComplexOrRealType>
class KetForTargetingExpression {

	using PairType = std::pair<SizeType, ComplexOrRealType>;
	using SumStruct = std::vector<PairType>;

public:

	using FactorForTargetingExpressionType = FactorForTargetingExpression<ComplexOrRealType>;

	enum class Kind {U, X, P, M, S};

	KetForTargetingExpression() : factor_(1.0), kind_(Kind::U) {}

	KetForTargetingExpression(const std::string& str) : factor_(1.0), kind_(Kind::U)
	{
		set(str);
	}

	void set(const std::string& str)
	{
		str_ = str;
		setKind();
	}

	void setFactor(const ComplexOrRealType& val)
	{
		factor_.set(val);
	}

	void multiply(const std::string& op)
	{
		if (kind_ == Kind::M) {
			err("KetForTargetingExpression: multiply one op at a time\n");
		}

		op_ = op;
		kind_ = Kind::M;
	}

	void multiply(const ComplexOrRealType& val)
	{
		factor_.multiply(val);
	}

	void sum(const KetForTargetingExpression& other)
	{
		const ComplexOrRealType& factor = other.factor_.value();
		if (kind_ == Kind::S) {
			assert(pVectors_to_sum_.size() > 0);
			factor_.push(factor);
		} else {
			assert(pVectors_to_sum_.empty());
			factor_.set(factor);
			kind_ = Kind::S;
		}

		if (other.pIndex_ < 0) {
			err("sum of kets internal error pindex < 0\n");
		}

		pVectors_to_sum_.push_back(other.pIndex_);
	}

	// constant member below

	std::string toString(const std::vector<std::string>& ops) const
	{
		PsimagLite::String s;

		PsimagLite::String f = factor_.toString();

		const SizeType n = ops.size();
		for (SizeType i = 0; i < n; ++i)
			s += ops[i] + "*";
		return f + s + this->name();
	}

	std::string name() const
	{
		if (kind_ == Kind::M) {
			assert(!op_.empty());
			return  op_ + "*" + str_;
		} else if (kind_ == Kind::S) {
			SizeType n = pVectors_to_sum_.size();
			assert(n > 0);

			std::string p0PlusP1 = "|P" + ttos(pVectors_to_sum_[0]) + ">";
			for (SizeType i = 1; i < n; ++i) {
				p0PlusP1 += "+|P" + ttos(pVectors_to_sum_[i]) + ">";
			}

			return p0PlusP1;
		}

		return str_;
	}

	int pIndex() const
	{
		return pIndex_;

	}

	Kind kind() const {return kind_; }

	const std::string op() const { return op_; }

	bool isSummable() const
	{
		return (kind_ == Kind::P);
	}

	bool canSumBeFinished() const
	{
		return (kind_ == Kind::S && pVectors_to_sum_.size() > 1);
	}

        ComplexOrRealType factor() const
	{
		return factor_.value();
	}

	SumStruct fillSumStruct() const
	{
		assert(pVectors_to_sum_.size() == 2);
		SizeType ind0 = pVectors_to_sum_[0];
		SizeType ind1 = pVectors_to_sum_[1];

		assert(ind0 != ind1);
		FactorForTargetingExpressionType f = factor_;
		if (ind0 > ind1) {
			std::swap(ind0, ind1);
			f.swap();
		}

		PairType p0(ind0, f.value(0));
		PairType p1(ind1, f.value(1));
		return {p0, p1};
	}
private:

	void setKind()
	{
		SizeType last = str_.length();
		if (last < 4) return;

		if (str_ == "|gs>") {
			kind_ = Kind::X;
			return;
		}

		--last;
		if (str_.substr(0, 2) == "|P" && str_[last] == '>') {
			pIndex_ = PsimagLite::atoi(str_.substr(2, last - 2));
			kind_ = Kind::P;
		} else {
			kind_ = Kind::U;
		}
	}

	std::string str_;
	FactorForTargetingExpressionType factor_;
	Kind kind_;
	int pIndex_ = -1;
	std::string op_;
	std::vector<SizeType> pVectors_to_sum_;
};
}
#endif // KETFORTARGETINGEXPRESSION_HH
