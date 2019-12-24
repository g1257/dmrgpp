#ifndef HAMILTONIANABSTRACT_H
#define HAMILTONIANABSTRACT_H

#include "Vector.h"
#include "ProgramGlobals.h"

namespace Dmrg {

template<typename SuperGeometryType>
class HamiltonianAbstract {

	typedef PsimagLite::Vector<SizeType>::Type VectorSizeType;
	typedef PsimagLite::Vector<VectorSizeType>::Type VectorVectorSizeType;

public:

	HamiltonianAbstract(const SuperGeometryType& superGeometry,
	                    SizeType smax,
	                    SizeType emin,
	                    const VectorSizeType& block)
	    : block_(block)
	{
		VectorSizeType v(2, 0);
		SizeType n = block.size();
		for (SizeType i = 0; i < n; ++i) {
			for (SizeType j = i + 1; j < n; ++j) {
				v[0] = block[i];
				v[1] = block[j];
				if (!superGeometry.connected(smax, emin, v)) continue;

				ProgramGlobals::ConnectionEnum type = superGeometry.connectionKind(smax, v);

				if (type == ProgramGlobals::ConnectionEnum::SYSTEM_SYSTEM ||
				        type == ProgramGlobals::ConnectionEnum::ENVIRON_ENVIRON) continue;

				data_.push_back(v);
			}
		}

		superGeometry.addSuperConnections(data_, smax, emin);
	}

	SizeType items() const { return data_.size(); }

	VectorSizeType item(SizeType ind) const
	{
		assert(ind < data_.size());
		return data_[ind];
	}

	const VectorSizeType& block() const { return block_; }

private:

	const VectorSizeType& block_;
	VectorVectorSizeType data_;
};
}
#endif // HAMILTONIANABSTRACT_H
