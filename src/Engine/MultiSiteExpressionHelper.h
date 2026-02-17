#ifndef MULTISITEEXPRESSIONHELPER_H
#define MULTISITEEXPRESSIONHELPER_H
#include "DmrgSerializer.h"
#include "Vector.h"

namespace Dmrg {

template <typename LeftRightSuperType, typename VectorWithOffsetType_>
class MultiSiteExpressionHelper {

public:

	using VectorWithOffsetType = VectorWithOffsetType_;
	using DmrgSerializerType   = DmrgSerializer<LeftRightSuperType, VectorWithOffsetType>;
	typedef
	    typename PsimagLite::Vector<DmrgSerializerType const*>::Type VectorDmrgSerializerType;
	using VectorVectorWithOffsetType = typename PsimagLite::Vector<VectorWithOffsetType>::Type;
	using BasisWithOperatorsType     = typename LeftRightSuperType::BasisWithOperatorsType;
	using FermionSignType            = typename DmrgSerializerType::FermionSignType;

	MultiSiteExpressionHelper(SizeType n)
	    : vds_(n, nullptr)
	{ }

	~MultiSiteExpressionHelper()
	{
		const SizeType n = vds_.size();
		for (SizeType i = 0; i < n; ++i) {
			delete vds_[i];
			vds_[i] = nullptr;
		}
	}

	void push(DmrgSerializerType const* ds, const VectorWithOffsetType& psi)
	{
		SizeType coo = ds->centerOfOrthogonality();
		assert(coo > 0);
		--coo;

		delete vds_[coo];
		vds_[coo] = ds;
		vgs_[coo] = psi;
	}

private:

	VectorDmrgSerializerType   vds_;
	VectorVectorWithOffsetType vgs_;
};
}
#endif // MULTISITEEXPRESSIONHELPER_H
