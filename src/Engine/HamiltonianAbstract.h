#ifndef HAMILTONIANABSTRACT_H
#define HAMILTONIANABSTRACT_H

#include "OneLink.hh"
#include "ProgramGlobals.h"
#include "Vector.h"

namespace Dmrg {

template <typename SuperGeometryType> class HamiltonianAbstract {
	using VectorSizeType = std::vector<SizeType>;
	using VectorVectorSizeType = std::vector<VectorSizeType>;
	using ComplexOrRealType = typename SuperGeometryType::ComplexOrRealType;
	using RealType = typename PsimagLite::Real<ComplexOrRealType>::Type;
	using OneLinkType = OneLink<ComplexOrRealType>;

public:

	HamiltonianAbstract(const SuperGeometryType& superGeometry,
	                    SizeType smax,
	                    SizeType emin,
	                    const VectorSizeType& block)
	    : superGeometry_(superGeometry)
	    , smax_(smax)
	    , emin_(emin)
	    , block_(block)
	{
		VectorSizeType v(2, 0);
		SizeType n = block.size();
		for (SizeType i = 0; i < n; ++i) {
			for (SizeType j = i + 1; j < n; ++j) {
				v[0] = block[i];
				v[1] = block[j];
				if (!superGeometry.connected(smax, emin, v))
					continue;

				ProgramGlobals::ConnectionEnum type
				    = superGeometry.connectionKind(smax, v);

				if (type == ProgramGlobals::ConnectionEnum::SYSTEM_SYSTEM
				    || type == ProgramGlobals::ConnectionEnum::ENVIRON_ENVIRON)
					continue;

				data_.push_back(v);
			}
		}

		superGeometry.addSuperConnections(data_, smax, emin);
	}

	ComplexOrRealType connectionValue(const VectorSizeType& hItems,
	                                  const OneLinkType& oneLink,
	                                  SizeType termIndexForGeom,
	                                  const RealType& targetTime)
	{
		ComplexOrRealType value
		    = superGeometry_(smax_, emin_, hItems, oneLink.orbs, termIndexForGeom);
		SizeType site = findSite(hItems);
		oneLink.modifier(value, targetTime, site);
		return value;
	}

	SizeType items() const { return data_.size(); }

	VectorSizeType item(SizeType ind) const
	{
		assert(ind < data_.size());
		return data_[ind];
	}

	const VectorSizeType& block() const { return block_; }

	const SuperGeometryType& superGeometry() const { return superGeometry_; }

private:

	SizeType findSite(const VectorSizeType& hItems) const
	{
		ProgramGlobals::ConnectionEnum type = superGeometry_.connectionKind(smax_, hItems);
		assert(type == ProgramGlobals::ConnectionEnum::SYSTEM_ENVIRON
		       || type == ProgramGlobals::ConnectionEnum::ENVIRON_SYSTEM);
		return (type == ProgramGlobals::ConnectionEnum::SYSTEM_ENVIRON) ? hItems[0]
		                                                                : hItems[1];
	}

	const SuperGeometryType& superGeometry_;
	const SizeType smax_;
	const SizeType emin_;
	const VectorSizeType& block_;
	VectorVectorSizeType data_;
};
}
#endif // HAMILTONIANABSTRACT_H
