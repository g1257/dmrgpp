#ifndef SUPERGEOMETRY_H
#define SUPERGEOMETRY_H
#include "Geometry/Geometry.h"
#include "ProgramGlobals.h"

namespace Dmrg {

template<typename GeometryType>
class SuperGeometry {

	typedef typename GeometryType::ComplexOrRealType ComplexOrRealType;
	typedef typename GeometryType::VectorSizeType VectorSizeType;
	typedef typename GeometryType::AdditionalDataType AdditionalDataType;

public:

	SuperGeometry(const GeometryType& geometry)
	    : geometry_(geometry)
	{}

	const GeometryType& geometry() const { return geometry_; }

	ComplexOrRealType operator()(SizeType smax,
	                             SizeType emin,
	                             const VectorSizeType& hItems,
	                             const VectorSizeType& edofs,
	                             SizeType term) const
	{
		checkVectorHasTwoEntries(hItems);
		checkVectorHasTwoEntries(edofs);
		return geometry_(smax, emin, hItems[0], edofs[0], hItems[1], edofs[1], term);
	}

	bool connected(SizeType smax,SizeType emin, const VectorSizeType& hItems) const
	{
		checkVectorHasTwoEntries(hItems);

		return geometry_.connected(smax, emin, hItems[0], hItems[1]);
	}

	typename ProgramGlobals::ConnectionEnum connectionKind(SizeType smax,
	                                                       const VectorSizeType& hItems) const
	{
		checkVectorHasTwoEntries(hItems);

		return geometry_.connectionKind(smax, hItems[0], hItems[1]);
	}

	void fillAdditionalData(AdditionalDataType& additionalData,
	                        SizeType term,
	                        const VectorSizeType& hItems) const
	{
		checkVectorHasTwoEntries(hItems);
		return geometry_.fillAdditionalData(additionalData, term, hItems[0], hItems[1]);
	}

private:

	static void checkVectorHasTwoEntries(const VectorSizeType& hItems)
	{
		if (hItems.size() != 2)
			err("SuperGeometry unimplemented\n");
	}

	const GeometryType& geometry_;
};
}
#endif // SUPERGEOMETRY_H
