#ifndef GEOMETRY_EX_H
#define GEOMETRY_EX_H

#include "Vector.h"
#include "String.h"

namespace PsimagLite {

template<typename GeometryTermType>
class GeometryExBase {

	typedef typename GeometryTermType::RealType RealType;

public:

	typedef typename Vector<RealType>::Type VectorRealType;

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
};

}

#ifndef USE_MS_GEOMETRY
namespace PsimagLite {
template<typename GeometryTermType>
class GeometryEx  : public GeometryExBase<GeometryTermType> {

	typedef typename Vector<GeometryTermType*>::Type VectorGeometryTermType;

public:

	GeometryEx(SizeType linSize,VectorGeometryTermType& terms)
	{}
};

}

#else

#include "msGeometry.h"

namespace PsimagLite {

template<typename GeometryTermType>
class GeometryEx  : public GeometryExBase<GeometryTermType> {

	typedef typename GeometryTermType::RealType RealType;
	typedef typename Vector<GeometryTermType*>::Type VectorGeometryTermType;

public:

	typedef typename Vector<RealType>::Type VectorRealType;

	GeometryEx(SizeType linSize,VectorGeometryTermType& terms)
	: linSize_(linSize),
	  terms_(terms)
	{
		SizeType ly = 2;
		SizeType lx = 2;
		SizeType meshPoints = 4;
		msGeometryInit2D(lx,0,0,ly,meshPoints);
	}

	~GeometryEx()
	{
		msGeometryEnd();
	}

	SizeType dimension(SizeType i = 0) const
	{
		assert(terms_.size() > i);
		String s = terms_[i]->label();

		if (s == "chain") return 1;

		if (s=="ladder" ||  s=="ladderbath") return 2;

		throw RuntimeError("GeometryEx: dimension\n");
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

private:

	SizeType linSize_;
	VectorGeometryTermType& terms_;
};

} // namespace PsimagLite
#endif

#endif // GEOMETRY_EX_H

