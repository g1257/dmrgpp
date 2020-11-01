#ifndef SDHSREINTERPRET_H
#define SDHSREINTERPRET_H
#include "Vector.h"

namespace Dmrg {

template<typename BraketType>
class SdhsReinterpret {

public:

	typedef PsimagLite::Vector<SizeType>::Type VectorSizeType;
	typedef typename BraketType::ModelType ModelType;
	typedef typename ModelType::OpaqueOp OpaqueOpType;
	typedef typename BraketType::RealType RealType;
	typedef typename BraketType::VectorAlgebraType VectorOperatorType;
	typedef typename BraketType::AlgebraType OperatorType;

	SdhsReinterpret(const BraketType& braket, const VectorSizeType& sites)
	   : forbidden_(false) // : braket_(braket), sites_(sites)
	{
		const SizeType n = sites.size();
		if (n == 0)
			err("SdhsReinterpret: 0 sites not allowed (FATAL)\n");
		if (n != braket.points())
			err("SdhsReinterpret: braket.points != sites.size()\n");

		for (SizeType i = 0; i < n; ++i) {
			PsimagLite::String opName = braket.opName(i);
			SizeType kind = braket.model().siteToAtomKind(sites[i]);
			OpaqueOpType opaque(opName);
			if (opaque.kindOfSite != kind) {
				forbidden_ = true;
				break;
			}

			ops_.push_back(braket.op(i));
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
//	const BraketType& braket_;
//	const VectorSizeType& sites_;
};
}
#endif // SDHSREINTERPRET_H
