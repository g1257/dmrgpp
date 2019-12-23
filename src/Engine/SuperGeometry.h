#ifndef SUPERGEOMETRY_H
#define SUPERGEOMETRY_H
#include "Geometry/Geometry.h"
#include "ProgramGlobals.h"

namespace Dmrg {

template<typename GeometryType, typename SuperConnectorType>
class SuperGeometry {

	typedef typename GeometryType::ComplexOrRealType ComplexOrRealType;
	typedef typename GeometryType::VectorSizeType VectorSizeType;
	typedef typename GeometryType::AdditionalDataType AdditionalDataType;
	typedef typename PsimagLite::Vector<VectorSizeType>::Type VectorVectorSizeType;

public:

	SuperGeometry(const GeometryType& geometry, const SuperConnectorType& superc)
	    : geometry_(geometry), superc_(superc)
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

	SizeType overSize(SizeType blockSize) const
	{
		return blockSize*(blockSize/2 + 1); // + superc_.size();
	}

//	SizeType addSuperConnections(VectorVectorSizeType& data,
//	                             SizeType smax,
//	                             SizeType emin,
//	                             const VectorSizeType& block,
//	                             SizeType counter) const
//	{
//		SizeType c = superc_.addSuperConnections(data, smax, emin, block);
//		return c + counter;
//	}

private:

	static void checkVectorHasTwoEntries(const VectorSizeType& hItems)
	{
		if (hItems.size() != 2)
			err("SuperGeometry unimplemented\n");
	}

	const GeometryType& geometry_;
	const SuperConnectorType& superc_;
};
}
#endif // SUPERGEOMETRY_H
