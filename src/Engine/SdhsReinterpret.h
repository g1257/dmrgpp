#ifndef SDHSREINTERPRET_H
#define SDHSREINTERPRET_H
#include "Vector.h"

namespace Dmrg {

template <typename BraketType> class SdhsReinterpret {

public:

	typedef PsimagLite::Vector<SizeType>::Type VectorSizeType;
	typedef typename BraketType::ModelType ModelType;
	typedef typename ModelType::OpaqueOp OpaqueOpType;
	typedef typename BraketType::RealType RealType;
	typedef typename BraketType::VectorAlgebraType VectorOperatorType;
	typedef typename BraketType::AlgebraType OperatorType;
	typedef typename BraketType::OneOperatorSpecType OneOperatorSpecType;
	typedef typename OneOperatorSpecType::SiteSplit SiteSplitType;

	SdhsReinterpret(const BraketType& braket, const VectorSizeType& sites)
	    : forbidden_(false)
	{
		const SizeType n = sites.size();
		if (n == 0)
			err("SdhsReinterpret: 0 sites not allowed (FATAL)\n");
		if (n != braket.points())
			err("SdhsReinterpret: braket.points != sites.size()\n");

		PsimagLite::String str("<");
		str += braket.bra().toString();
		str += "|";
		for (SizeType i = 0; i < n; ++i) {
			PsimagLite::String opName = braket.opName(i);
			str += opName;
			SiteSplitType siteSplit = OneOperatorSpecType::extractSiteIfAny(opName);
			if (!siteSplit.hasSiteString)
				str += "[" + ttos(sites[i]) + "]";
			if (i < n - 1)
				str += ";";
		}

		str += "|" + braket.ket().toString() + ">";
		BraketType braket2(braket.model(), str);
		for (SizeType i = 0; i < n; ++i) {
			ops_.push_back(braket2.op(i));
			if (ops_[i].isEmpty())
				forbidden_ = true;
		}
	}

	bool forbidden() const { return forbidden_; }

	static RealType forbiddenValue() { return -100; }

	const OperatorType& op(SizeType ind) const
	{
		assert(ind < ops_.size());
		return ops_[ind];
	}

private:

	bool forbidden_;
	VectorOperatorType ops_;
};
}
#endif // SDHSREINTERPRET_H
