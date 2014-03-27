#ifndef GEOMETRY_EX_H
#define GEOMETRY_EX_H

#include "Vector.h"
#include "String.h"

#ifndef USE_MS_GEOMETRY

namespace PsimagLite {

template<typename RealType>
class GeometryEx {

public:

	typedef typename Vector<RealType>::Type VectorRealType;

	GeometryEx(SizeType) {}

	SizeType dimension(SizeType i = 0) const
	{
		throw RuntimeError("GeometryEx: dimension\n");
	}

	void index2Rvector(SizeType r,VectorRealType& rvector) const
	{
		throw RuntimeError("GeometryEx: index2Rvector\n");
	}

	void index2Kvector(SizeType k,VectorRealType& kvector) const
	{
		throw RuntimeError("GeometryEx: index2Kvector\n");
	}

	//! Number of symmetry operations for this K Geometry.
	SizeType nGroupK() const
	{
		throw RuntimeError("GeometryEx: nGroupK\n");
	}

	SizeType ickequ(SizeType j,SizeType op) const
	{
		throw RuntimeError("GeometryEx: ickequ\n");
	}

	void getMeshVector(VectorRealType& kvector,SizeType k) const
	{
		throw RuntimeError("GeometryEx: getMeshVector\n");
	}

	SizeType sizeOfMesh() const
	{
		throw RuntimeError("GeometryEx: sizeOfMesh\n");
	}
};

}

#else

#include "msGeometry.h"

namespace PsimagLite {

template<typename RealType>
class GeometryEx  {

public:

	typedef typename Vector<RealType>::Type VectorRealType;

	GeometryEx(SizeType meshPoints)
	{
		SizeType ly = 2;
		SizeType lx = 2;
		std::cerr<<"WARNING: GeometryEx(): lattice of 2x2 hard wired\n";
		msGeometryInit2D(lx,0,0,ly,meshPoints);
	}

	~GeometryEx()
	{
		msGeometryEnd();
	}

	SizeType dimension(SizeType i = 0) const
	{
		return msGeometryDim();
	}

	void index2Rvector(SizeType r,VectorRealType& rvector) const
	{
		return msGeometryIndex2Rvector(r,rvector);
	}

	void index2Kvector(SizeType k,VectorRealType& kvector) const
	{
		return msGeometryIndex2Kvector(k,kvector);
	}

	//! Number of symmetry operations for this K Geometry.
	SizeType nGroupK() const
	{
		return msGeometryNgroupK();
	}

	SizeType ickequ(SizeType j,SizeType op) const
	{
		return msGeometryIckequ(j,op);
	}

	void getMeshVector(VectorRealType& kvector,SizeType k) const
	{
		return msGeometryGetMeshVector(kvector,k);
	}

	SizeType sizeOfMesh() const
	{
		return msGeometrySizeOfMesh();
	}
};

} // namespace PsimagLite
#endif

#endif // GEOMETRY_EX_H

