#ifndef MULTISITEEXPRESSIONHELPER_H
#define MULTISITEEXPRESSIONHELPER_H
#include "DmrgSerializer.h"
#include "Vector.h"

namespace Dmrg {

template<typename LeftRightSuperType, typename VectorWithOffsetType_>
class MultiSiteExpressionHelper {

public:

	typedef VectorWithOffsetType_ VectorWithOffsetType;
	typedef DmrgSerializer<LeftRightSuperType, VectorWithOffsetType> DmrgSerializerType;
	typedef typename PsimagLite::Vector<DmrgSerializerType const*>::Type VectorDmrgSerializerType;
	typedef typename PsimagLite::Vector<VectorWithOffsetType>::Type VectorVectorWithOffsetType;
	typedef typename LeftRightSuperType::BasisWithOperatorsType BasisWithOperatorsType;
	typedef typename DmrgSerializerType::FermionSignType FermionSignType;

	MultiSiteExpressionHelper(SizeType n) : vds_(n, nullptr) {}

	~MultiSiteExpressionHelper()
	{
		const SizeType n = vds_.size();
		for (SizeType i = 0; i < n; ++i) {
			delete vds_[i];
			vds_[i] = nullptr;
		}
	}

	void push(DmrgSerializerType const* ds,
	          const VectorWithOffsetType& psi)
	{
		SizeType coo = ds->centerOfOrthogonality();
		assert(coo > 0);
		--coo;

		delete vds_[coo];
		vds_[coo] = ds;
		vgs_[coo] = psi;
	}

private:

	VectorDmrgSerializerType vds_;
	VectorVectorWithOffsetType vgs_;
};
}
#endif // MULTISITEEXPRESSIONHELPER_H
